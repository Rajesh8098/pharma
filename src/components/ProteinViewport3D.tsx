import React, { useEffect, useRef } from 'react';
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';
import {
  AtomRecord,
  ColorSchemeMode,
  HydrogenBond,
  ProteinStructureModel,
  RepresentationMode,
  ResidueRecord,
} from '../types/protein';
import { getPlddtColorHex, HYDROPHOBICITY_SCALE } from '../utils/proteinParser';

export type CameraPresetType = 'overview' | 'top' | 'side' | 'focus-selected' | 'focus-hbonds';

interface ProteinViewport3DProps {
  model: ProteinStructureModel;
  representation: RepresentationMode;
  colorScheme: ColorSchemeMode;
  showHydrogenBonds: boolean;
  filteredHydrogenBonds: HydrogenBond[];
  selectedResidueKey: string | null;
  hoveredResidueKey: string | null;
  activeHBondId: string | null;
  lightBackground: boolean;
  explodedScale: number;
  cameraPresetTrigger: { preset: CameraPresetType; seq: number };
  onSelectResidue: (residueKey: string) => void;
  onHoverResidue: (residueKey: string | null, atomInfo?: AtomRecord | null) => void;
}

// Pre-allocated module-scope scratch objects (Zero GC inside render/build loops)
const _vecA = new THREE.Vector3();
const _vecB = new THREE.Vector3();
const _mid = new THREE.Vector3();
const _dir = new THREE.Vector3();
const _up = new THREE.Vector3(0, 1, 0);
const _quat = new THREE.Quaternion();
const _scale = new THREE.Vector3();
const _mat4 = new THREE.Matrix4();
const _color = new THREE.Color();

const ELEMENT_COLORS: Record<string, number> = {
  C: 0x94a3b8,
  N: 0x3b82f6,
  O: 0xef4444,
  S: 0xeab308,
  FE: 0xf97316,
  ZN: 0xa855f7,
  MG: 0x22c55e,
  P: 0xf59e0b,
};

const CHAIN_PALETTE = [0x38bdf8, 0x34d399, 0xf472b6, 0xfbbf24, 0xa78bfa];

const HBOND_CATEGORY_COLORS: Record<string, number> = {
  'backbone-helix': 0x38bdf8, // Cyan
  'backbone-sheet': 0xc084fc, // Violet
  sidechain: 0x10b981,        // Emerald
  ligand: 0xf59e0b,           // Amber
};

function getResidueColorHexNumber(
  res: ResidueRecord,
  colorScheme: ColorSchemeMode,
  chainIndex: number
): number {
  if (res.isLigand) return 0xf59e0b;

  switch (colorScheme) {
    case 'plddt-bfactor':
      return new THREE.Color(getPlddtColorHex(res.plddtOrBfactor)).getHex();
    case 'secondary-structure':
      if (res.secondaryStructure === 'HELIX') return 0x06b6d4; // Teal/Cyan Helix
      if (res.secondaryStructure === 'SHEET') return 0xf59e0b; // Gold Beta Sheet
      return 0x94a3b8;                                         // Slate Coil
    case 'chain':
      return CHAIN_PALETTE[chainIndex % CHAIN_PALETTE.length];
    case 'hydrophobicity': {
      const hydro = HYDROPHOBICITY_SCALE[res.resName] ?? 0;
      // Map [-4.5, +4.5] from hydrophilic cobalt blue (0x2563eb) to hydrophobic amber-orange (0xea580c)
      const t = Math.max(0, Math.min(1, (hydro + 4.5) / 9.0));
      const c1 = new THREE.Color(0x2563eb);
      const c2 = new THREE.Color(0xea580c);
      return c1.lerp(c2, t).getHex();
    }
    case 'element':
    default:
      return 0x64748b;
  }
}

/**
 * Generates a procedural radial diffraction / goniometer dial texture and bump map (Zero external URLs)
 */
function createProceduralDiffractionTexture(maxAnisotropy: number): {
  colorMap: THREE.CanvasTexture;
  bumpMap: THREE.CanvasTexture;
} {
  const size = 512;
  const canvas = document.createElement('canvas');
  canvas.width = size;
  canvas.height = size;
  const ctx = canvas.getContext('2d')!;

  ctx.fillStyle = '#0f172a';
  ctx.fillRect(0, 0, size, size);

  const cx = size / 2;
  const cy = size / 2;

  // Concentric Bragg diffraction resolution rings
  const radii = [60, 115, 165, 210, 244];
  radii.forEach((r, idx) => {
    ctx.beginPath();
    ctx.arc(cx, cy, r, 0, Math.PI * 2);
    ctx.strokeStyle = idx % 2 === 0 ? 'rgba(56, 189, 248, 0.35)' : 'rgba(148, 163, 184, 0.22)';
    ctx.lineWidth = idx === radii.length - 1 ? 3 : 1.5;
    ctx.stroke();
  });

  // 72 radial degree ticks around the perimeter
  for (let deg = 0; deg < 360; deg += 5) {
    const rad = (deg * Math.PI) / 180;
    const inner = deg % 30 === 0 ? 192 : deg % 15 === 0 ? 204 : 214;
    const outer = 224;
    ctx.beginPath();
    ctx.moveTo(cx + Math.cos(rad) * inner, cy + Math.sin(rad) * inner);
    ctx.lineTo(cx + Math.cos(rad) * outer, cy + Math.sin(rad) * outer);
    ctx.strokeStyle = deg % 30 === 0 ? 'rgba(56, 189, 248, 0.65)' : 'rgba(148, 163, 184, 0.28)';
    ctx.lineWidth = deg % 30 === 0 ? 2 : 1;
    ctx.stroke();
  }

  const colorMap = new THREE.CanvasTexture(canvas);
  colorMap.colorSpace = THREE.SRGBColorSpace;
  colorMap.anisotropy = maxAnisotropy;

  // Bump map canvas
  const bumpCanvas = document.createElement('canvas');
  bumpCanvas.width = size;
  bumpCanvas.height = size;
  const bCtx = bumpCanvas.getContext('2d')!;
  bCtx.fillStyle = '#808080';
  bCtx.fillRect(0, 0, size, size);
  radii.forEach((r) => {
    bCtx.beginPath();
    bCtx.arc(cx, cy, r, 0, Math.PI * 2);
    bCtx.strokeStyle = '#ffffff';
    bCtx.lineWidth = 2;
    bCtx.stroke();
  });
  const bumpMap = new THREE.CanvasTexture(bumpCanvas);
  bumpMap.colorSpace = THREE.NoColorSpace;
  bumpMap.anisotropy = maxAnisotropy;

  return { colorMap, bumpMap };
}

/**
 * Generates a soft radial contact shadow texture so the structure & pedestal sit grounded
 */
function createRadialShadowTexture(): THREE.CanvasTexture {
  const size = 256;
  const canvas = document.createElement('canvas');
  canvas.width = size;
  canvas.height = size;
  const ctx = canvas.getContext('2d')!;
  const grad = ctx.createRadialGradient(size / 2, size / 2, 10, size / 2, size / 2, size / 2);
  grad.addColorStop(0, 'rgba(0, 0, 0, 0.68)');
  grad.addColorStop(0.55, 'rgba(0, 0, 0, 0.25)');
  grad.addColorStop(1, 'rgba(0, 0, 0, 0.0)');
  ctx.fillStyle = grad;
  ctx.fillRect(0, 0, size, size);
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.SRGBColorSpace;
  return tex;
}

/**
 * Builds the 7-part Cryo-EM / X-Ray Diffraction Goniometer Pedestal Assembly (Composite THREE.Group)
 * Pivot is at bottom-center flush with the floor plane.
 */
function createGoniometerStageGroup(maxAnisotropy: number): THREE.Group {
  const group = new THREE.Group();
  group.name = 'CryoEMGoniometerStage';

  const { colorMap, bumpMap } = createProceduralDiffractionTexture(maxAnisotropy);
  const shadowTex = createRadialShadowTexture();

  // 1. Soft radial contact shadow plane at y = 0.015
  const shadowGeo = new THREE.CircleGeometry(21, 48);
  shadowGeo.rotateX(-Math.PI / 2);
  const shadowMat = new THREE.MeshBasicMaterial({
    map: shadowTex,
    transparent: true,
    opacity: 0.48,
    depthWrite: false,
  });
  const shadowMesh = new THREE.Mesh(shadowGeo, shadowMat);
  shadowMesh.position.y = 0.015;
  group.add(shadowMesh);

  // 2. Primary Silhouette: Outer Beveled Titanium Goniometer Ring via LatheGeometry
  const lathePoints = [
    new THREE.Vector2(13.2, 0.0),
    new THREE.Vector2(15.6, 0.0),
    new THREE.Vector2(16.0, 0.35),
    new THREE.Vector2(15.5, 0.85),
    new THREE.Vector2(13.6, 0.85),
    new THREE.Vector2(13.2, 0.5),
    new THREE.Vector2(13.2, 0.0),
  ];
  const outerRingGeo = new THREE.LatheGeometry(lathePoints, 64);
  const outerRingMat = new THREE.MeshPhysicalMaterial({
    color: 0x334155,
    roughness: 0.28,
    metalness: 0.86,
    clearcoat: 0.3,
    bumpMap,
    bumpScale: 0.03,
  });
  const outerRingMesh = new THREE.Mesh(outerRingGeo, outerRingMat);
  outerRingMesh.receiveShadow = true;
  group.add(outerRingMesh);

  // 3. Secondary Read: Calibrated Diffraction Plate Disc with procedural Bragg markings
  const plateGeo = new THREE.CylinderGeometry(13.2, 13.2, 0.24, 64);
  const plateMat = new THREE.MeshStandardMaterial({
    color: 0xffffff,
    map: colorMap,
    bumpMap,
    bumpScale: 0.015,
    roughness: 0.42,
    metalness: 0.22,
  });
  const plateMesh = new THREE.Mesh(plateGeo, plateMat);
  plateMesh.position.y = 0.22;
  plateMesh.receiveShadow = true;
  group.add(plateMesh);

  // 4. Secondary Read: Inner Anodized Vernier Azimuth Ring
  const vernierGeo = new THREE.TorusGeometry(9.4, 0.14, 16, 64);
  vernierGeo.rotateX(Math.PI / 2);
  const vernierMat = new THREE.MeshStandardMaterial({
    color: 0x0284c7,
    roughness: 0.22,
    metalness: 0.78,
  });
  const vernierMesh = new THREE.Mesh(vernierGeo, vernierMat);
  vernierMesh.position.y = 0.38;
  group.add(vernierMesh);

  // 5. Accent Details: 4 Cardinal Beveled Alignment Pylons (ExtrudeGeometry) + Brass Index Pins
  const pylonShape = new THREE.Shape();
  pylonShape.moveTo(-0.45, 0);
  pylonShape.lineTo(0.45, 0);
  pylonShape.lineTo(0.32, 1.15);
  pylonShape.lineTo(-0.32, 1.15);
  pylonShape.closePath();

  const pylonGeo = new THREE.ExtrudeGeometry(pylonShape, {
    steps: 1,
    depth: 0.9,
    bevelEnabled: true,
    bevelThickness: 0.08,
    bevelSize: 0.06,
    bevelSegments: 2,
  });
  pylonGeo.center();

  const pylonMat = new THREE.MeshStandardMaterial({
    color: 0x475569,
    roughness: 0.62,
    metalness: 0.45,
  });

  const pinGeo = new THREE.CylinderGeometry(0.14, 0.14, 0.4, 16);
  const pinMat = new THREE.MeshStandardMaterial({
    color: 0xf59e0b,
    roughness: 0.18,
    metalness: 0.92,
  });

  for (let k = 0; k < 4; k++) {
    const angle = (k * Math.PI) / 2;
    const pylon = new THREE.Mesh(pylonGeo, pylonMat);
    pylon.position.set(Math.cos(angle) * 14.5, 0.75, Math.sin(angle) * 14.5);
    pylon.rotation.y = -angle;
    group.add(pylon);

    const pin = new THREE.Mesh(pinGeo, pinMat);
    pin.position.set(Math.cos(angle) * 14.5, 1.4, Math.sin(angle) * 14.5);
    group.add(pin);
  }

  return group;
}

export const ProteinViewport3D: React.FC<ProteinViewport3DProps> = ({
  model,
  representation,
  colorScheme,
  showHydrogenBonds,
  filteredHydrogenBonds,
  selectedResidueKey,
  hoveredResidueKey,
  activeHBondId,
  lightBackground,
  explodedScale,
  cameraPresetTrigger,
  onSelectResidue,
  onHoverResidue,
}) => {
  const canvasMountRef = useRef<HTMLDivElement | null>(null);

  // Store latest callbacks in refs so React updates never re-initialize the WebGL context
  const onSelectRef = useRef(onSelectResidue);
  const onHoverRef = useRef(onHoverResidue);
  onSelectRef.current = onSelectResidue;
  onHoverRef.current = onHoverResidue;

  // Store Three.js scene handles in ref
  const threeCtxRef = useRef<{
    renderer: THREE.WebGLRenderer;
    scene: THREE.Scene;
    camera: THREE.PerspectiveCamera;
    controls: OrbitControls;
    proteinGroup: THREE.Group;
    hbondGroup: THREE.Group;
    beaconGroup: THREE.Group;
    stageGroup: THREE.Group;
    atomInstancedMesh: THREE.InstancedMesh | null;
    instanceToAtomMap: AtomRecord[];
    centroid: THREE.Vector3;
    boundingRadius: number;
    desiredCamPos: THREE.Vector3 | null;
    desiredCamTarget: THREE.Vector3 | null;
    requestRender: () => void;
    isSoftwareGl: boolean;
  } | null>(null);

  // 1. Initialize WebGLRenderer, Studio Environment, Lights, and Calibration Stage ONCE
  useEffect(() => {
    const container = canvasMountRef.current;
    if (!container) return;

    const isSoftwareGl = (() => {
      const gl = document.createElement('canvas').getContext('webgl2');
      if (!gl) return true;
      const ext = gl.getExtension('WEBGL_debug_renderer_info');
      const param = ext ? ext.UNMASKED_RENDERER_WEBGL : gl.RENDERER;
      const gpu = String(gl.getParameter(param));
      gl.getExtension('WEBGL_lose_context')?.loseContext();
      return /SwiftShader|llvmpipe|softpipe|Software|Basic Render/i.test(gpu);
    })();

    const width = container.clientWidth || window.innerWidth;
    const height = container.clientHeight || window.innerHeight;

    const renderer = new THREE.WebGLRenderer({
      antialias: !isSoftwareGl,
      powerPreference: 'high-performance',
      preserveDrawingBuffer: true,
    });
    renderer.outputColorSpace = THREE.SRGBColorSpace;
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = 1.12;
    renderer.setPixelRatio(isSoftwareGl ? 1 : Math.min(window.devicePixelRatio, 2));
    renderer.setSize(width, height);
    renderer.shadowMap.enabled = !isSoftwareGl;
    renderer.shadowMap.type = THREE.PCFSoftShadowMap;

    container.innerHTML = '';
    container.appendChild(renderer.domElement);

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0b0f19);
    scene.fog = new THREE.FogExp2(0x0b0f19, 0.0085);

    // Studio Environment Reflection (RoomEnvironment + PMREMGenerator)
    const pmremGenerator = new THREE.PMREMGenerator(renderer);
    scene.environment = pmremGenerator.fromScene(new RoomEnvironment(), 0.04).texture;
    pmremGenerator.dispose();

    const camera = new THREE.PerspectiveCamera(45, width / height, 0.5, 400);
    camera.position.set(0, 6, 38);

    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = !isSoftwareGl;
    controls.dampingFactor = 0.07;
    controls.minDistance = 6;
    controls.maxDistance = 120;
    controls.target.set(0, 0, 0);

    // 3-Point Studio Lighting Rig (Warm Key Light + Cool Hemisphere Fill + Cyan Rim Light)
    const hemiLight = new THREE.HemisphereLight(0xdbeafe, 0x1e293b, 0.85);
    hemiLight.position.set(0, 35, 0);
    scene.add(hemiLight);

    const keyLight = new THREE.DirectionalLight(0xfff7ed, 1.65);
    keyLight.position.set(18, 28, 22);
    keyLight.castShadow = !isSoftwareGl;
    keyLight.shadow.mapSize.set(1024, 1024);
    keyLight.shadow.camera.near = 2;
    keyLight.shadow.camera.far = 90;
    const d = 24;
    keyLight.shadow.camera.left = -d;
    keyLight.shadow.camera.right = d;
    keyLight.shadow.camera.top = d;
    keyLight.shadow.camera.bottom = -d;
    keyLight.shadow.bias = -0.0005;
    keyLight.shadow.normalBias = 0.02;
    scene.add(keyLight);

    const rimLight = new THREE.DirectionalLight(0x38bdf8, 0.75);
    rimLight.position.set(-22, 12, -18);
    scene.add(rimLight);

    // Scene Hierarchy Groups
    const stageGroup = createGoniometerStageGroup(renderer.capabilities.getMaxAnisotropy());
    StagePositionDefault(stageGroup);
    scene.add(stageGroup);

    const proteinGroup = new THREE.Group();
    proteinGroup.name = 'MacromolecularAssembly';
    scene.add(proteinGroup);

    const hbondGroup = new THREE.Group();
    hbondGroup.name = 'HydrogenBondNetwork';
    scene.add(hbondGroup);

    const beaconGroup = new THREE.Group();
    beaconGroup.name = 'SelectedResidueBeacon';
    scene.add(beaconGroup);

    // On-demand coalesced rendering with software-WebGL fenceSync guard
    let frameId = 0;
    let renderBusy = false;
    let renderDirty = false;
    let lastTime = performance.now();

    function afterDraw(cb: () => void) {
      const gl = renderer.getContext() as WebGL2RenderingContext;
      const sync = gl.fenceSync ? gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0) : null;
      if (gl.flush) gl.flush();
      if (!sync) {
        frameId = window.setTimeout(cb, 80);
        return;
      }
      const poll = () => {
        const state = gl.getSyncParameter(sync, gl.SYNC_STATUS);
        if (state === gl.UNSIGNALED) {
          frameId = window.setTimeout(poll, 5);
        } else {
          gl.deleteSync(sync);
          frameId = window.setTimeout(cb, 80);
        }
      };
      frameId = window.setTimeout(poll, 5);
    }

    function requestRender() {
      if (renderBusy) {
        renderDirty = true;
        return;
      }
      renderBusy = true;
      frameId = requestAnimationFrame((now) => {
        renderDirty = false;
        const dt = Math.min(0.1, Math.max(0.001, (now - lastTime) / 1000));
        lastTime = now;

        const ctx = threeCtxRef.current;
        let stillMoving = false;

        if (ctx && ctx.desiredCamPos && ctx.desiredCamTarget) {
          if (isSoftwareGl) {
            camera.position.copy(ctx.desiredCamPos);
            controls.target.copy(ctx.desiredCamTarget);
            ctx.desiredCamPos = null;
            ctx.desiredCamTarget = null;
          } else {
            const lambda = 9.0;
            camera.position.x = THREE.MathUtils.damp(camera.position.x, ctx.desiredCamPos.x, lambda, dt);
            camera.position.y = THREE.MathUtils.damp(camera.position.y, ctx.desiredCamPos.y, lambda, dt);
            camera.position.z = THREE.MathUtils.damp(camera.position.z, ctx.desiredCamPos.z, lambda, dt);

            controls.target.x = THREE.MathUtils.damp(controls.target.x, ctx.desiredCamTarget.x, lambda, dt);
            controls.target.y = THREE.MathUtils.damp(controls.target.y, ctx.desiredCamTarget.y, lambda, dt);
            controls.target.z = THREE.MathUtils.damp(controls.target.z, ctx.desiredCamTarget.z, lambda, dt);

            if (
              camera.position.distanceToSquared(ctx.desiredCamPos) < 0.02 &&
              controls.target.distanceToSquared(ctx.desiredCamTarget) < 0.02
            ) {
              ctx.desiredCamPos = null;
              ctx.desiredCamTarget = null;
            } else {
              stillMoving = true;
            }
          }
        }

        controls.update();
        renderer.render(scene, camera);

        const done = () => {
          renderBusy = false;
          if (renderDirty || stillMoving) {
            requestRender();
          }
        };
        if (isSoftwareGl) afterDraw(done);
        else done();
      });
    }

    controls.addEventListener('start', () => {
      if (threeCtxRef.current) {
        threeCtxRef.current.desiredCamPos = null;
        threeCtxRef.current.desiredCamTarget = null;
      }
    });
    controls.addEventListener('change', () => {
      requestRender();
    });

    // Raycasting for 3D atom/residue click & hover
    const raycaster = new THREE.Raycaster();
    raycaster.layers.set(1);
    const pointer = new THREE.Vector2();
    let pointerDownPos = { x: 0, y: 0 };
    let lastHoverTime = 0;

    const handlePointerDown = (e: PointerEvent) => {
      pointerDownPos = { x: e.clientX, y: e.clientY };
    };

    const handlePointerUp = (e: PointerEvent) => {
      const dx = e.clientX - pointerDownPos.x;
      const dy = e.clientY - pointerDownPos.y;
      if (dx * dx + dy * dy > 36) return; // Ignore orbit drags

      const ctx = threeCtxRef.current;
      if (!ctx || !ctx.atomInstancedMesh) return;

      const rect = renderer.domElement.getBoundingClientRect();
      pointer.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
      pointer.y = -((e.clientY - rect.top) / rect.height) * 2 + 1;

      raycaster.setFromCamera(pointer, camera);
      const hits = raycaster.intersectObject(ctx.atomInstancedMesh, false);
      if (hits.length > 0 && hits[0].instanceId !== undefined) {
        const atom = ctx.instanceToAtomMap[hits[0].instanceId];
        if (atom) {
          onSelectRef.current(`${atom.chainId}:${atom.resSeq}`);
        }
      }
    };

    const handlePointerMove = (e: PointerEvent) => {
      const now = performance.now();
      if (now - lastHoverTime < 45) return; // Throttle hover raycast to ~22fps
      lastHoverTime = now;

      const ctx = threeCtxRef.current;
      if (!ctx || !ctx.atomInstancedMesh) return;

      const rect = renderer.domElement.getBoundingClientRect();
      pointer.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
      pointer.y = -((e.clientY - rect.top) / rect.height) * 2 + 1;

      raycaster.setFromCamera(pointer, camera);
      const hits = raycaster.intersectObject(ctx.atomInstancedMesh, false);
      if (hits.length > 0 && hits[0].instanceId !== undefined) {
        const atom = ctx.instanceToAtomMap[hits[0].instanceId];
        if (atom) {
          renderer.domElement.style.cursor = 'pointer';
          onHoverRef.current(`${atom.chainId}:${atom.resSeq}`, atom);
          return;
        }
      }
      renderer.domElement.style.cursor = 'grab';
      onHoverRef.current(null, null);
    };

    renderer.domElement.addEventListener('pointerdown', handlePointerDown);
    renderer.domElement.addEventListener('pointerup', handlePointerUp);
    renderer.domElement.addEventListener('pointermove', handlePointerMove);

    const handleResize = () => {
      if (!container) return;
      const w = container.clientWidth || window.innerWidth;
      const h = container.clientHeight || window.innerHeight;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
      requestRender();
    };
    window.addEventListener('resize', handleResize);

    threeCtxRef.current = {
      renderer,
      scene,
      camera,
      controls,
      proteinGroup,
      hbondGroup,
      beaconGroup,
      stageGroup,
      atomInstancedMesh: null,
      instanceToAtomMap: [],
      centroid: new THREE.Vector3(),
      boundingRadius: 15,
      desiredCamPos: null,
      desiredCamTarget: null,
      requestRender,
      isSoftwareGl,
    };

    if (!isSoftwareGl) {
      renderer.compile(scene, camera);
    }
    requestRender();

    return () => {
      cancelAnimationFrame(frameId);
      clearTimeout(frameId);
      window.removeEventListener('resize', handleResize);
      renderer.domElement.removeEventListener('pointerdown', handlePointerDown);
      renderer.domElement.removeEventListener('pointerup', handlePointerUp);
      renderer.domElement.removeEventListener('pointermove', handlePointerMove);
      controls.dispose();
      renderer.dispose();
    };
  }, []);

  // 2. Update Background & Fog when lightBackground toggles
  useEffect(() => {
    const ctx = threeCtxRef.current;
    if (!ctx) return;
    const bgHex = lightBackground ? 0xf1f5f9 : 0x0b0f19;
    ctx.scene.background = new THREE.Color(bgHex);
    ctx.scene.fog = new THREE.FogExp2(bgHex, 0.0075);
    ctx.requestRender();
  }, [lightBackground]);

  // 3. Rebuild Macromolecular 3D Geometry when model, representation, colorScheme, or explodedScale changes
  useEffect(() => {
    const ctx = threeCtxRef.current;
    if (!ctx) return;

    disposeGroupChildren(ctx.proteinGroup);
    ctx.atomInstancedMesh = null;
    ctx.instanceToAtomMap = [];

    const atoms = model.atoms;
    const residues = model.residues;
    if (atoms.length === 0) {
      ctx.requestRender();
      return;
    }

    // Compute geometric centroid and bounding radius so the structure is centered at (0, 0, 0)
    let cx = 0;
    let cy = 0;
    let cz = 0;
    let minY = Infinity;
    for (const a of atoms) {
      cx += a.x;
      cy += a.y;
      cz += a.z;
    }
    cx /= atoms.length;
    cy /= atoms.length;
    cz /= atoms.length;
    ctx.centroid.set(cx, cy, cz);

    let maxR2 = 64;
    for (const a of atoms) {
      const dx = a.x - cx;
      const dy = a.y - cy;
      const dz = a.z - cz;
      const r2 = dx * dx + dy * dy + dz * dz;
      if (r2 > maxR2) maxR2 = r2;
      if (dy < minY) minY = dy;
    }
    ctx.boundingRadius = Math.sqrt(maxR2);

    // Position the Goniometer Calibration Stage flush underneath the lowest atom
    ctx.stageGroup.position.set(0, Math.min(-9.5, minY - 2.6), 0);

    // Compute per-residue domain/secondary-structure exploded offset when explodedScale > 0
    const resOffsetMap = new Map<string, THREE.Vector3>();
    for (const res of residues) {
      if (explodedScale <= 0.001) {
        resOffsetMap.set(res.key, new THREE.Vector3(0, 0, 0));
      } else {
        const vx = res.caAtom.x - cx;
        const vy = res.caAtom.y - cy;
        const vz = res.caAtom.z - cz;
        const len = Math.sqrt(vx * vx + vy * vy + vz * vz) || 1;
        const factor =
          (res.isLigand ? 5.5 : res.secondaryStructure === 'HELIX' ? 3.5 : res.secondaryStructure === 'SHEET' ? 2.5 : 1.5) *
          explodedScale;
        resOffsetMap.set(
          res.key,
          new THREE.Vector3((vx / len) * factor, (vy / len) * factor, (vz / len) * factor)
        );
      }
    }

    // Helper for transformed atom position
    const getAtomWorldPos = (atom: AtomRecord, target: THREE.Vector3): THREE.Vector3 => {
      const off = resOffsetMap.get(`${atom.chainId}:${atom.resSeq}`);
      target.set(
        atom.x - cx + (off ? off.x : 0),
        atom.y - cy + (off ? off.y : 0),
        atom.z - cz + (off ? off.z : 0)
      );
      return target;
    };

    // A. Secondary Structure Ribbons (Cartoon / Backbone Trace)
    if (representation === 'cartoon-sidechain' || representation === 'backbone-trace') {
      const polyResidues = residues.filter((r) => !r.isLigand);
      let segStart = 0;

      while (segStart < polyResidues.length) {
        const firstRes = polyResidues[segStart];
        const ssType = firstRes.secondaryStructure;
        const chainId = firstRes.chainId;
        let segEnd = segStart;

        while (
          segEnd + 1 < polyResidues.length &&
          polyResidues[segEnd + 1].chainId === chainId &&
          (representation === 'backbone-trace' ||
            polyResidues[segEnd + 1].secondaryStructure === ssType)
        ) {
          segEnd++;
        }

        // Extend 1 residue on each side for seamless spline continuity
        const splineStart = Math.max(0, segStart - 1);
        const splineEnd = Math.min(polyResidues.length - 1, segEnd + 1);
        const points: THREE.Vector3[] = [];
        const segResidues: ResidueRecord[] = [];

        for (let k = splineStart; k <= splineEnd; k++) {
          if (polyResidues[k].chainId !== chainId) continue;
          const pt = new THREE.Vector3();
          getAtomWorldPos(polyResidues[k].caAtom, pt);
          points.push(pt);
          segResidues.push(polyResidues[k]);
        }

        if (points.length >= 2) {
          const curve = new THREE.CatmullRomCurve3(points, false, 'centripetal', 0.5);
          const tubularSegments = Math.max(12, points.length * 6);
          const radius =
            representation === 'backbone-trace'
              ? 0.28
              : ssType === 'HELIX'
              ? 0.62
              : ssType === 'SHEET'
              ? 0.52
              : 0.24;

          const tubeGeo = new THREE.TubeGeometry(curve, tubularSegments, radius, 12, false);

          // Apply per-vertex colors along the spline matching the residue color scheme
          const posAttr = tubeGeo.attributes.position;
          const colors = new Float32Array(posAttr.count * 3);
          const chainIdx = Math.max(0, model.chains.indexOf(chainId));

          for (let v = 0; v < posAttr.count; v++) {
            const vx = posAttr.getX(v);
            const vy = posAttr.getY(v);
            const vz = posAttr.getZ(v);

            // Find closest residue in this segment
            let bestRes = segResidues[0];
            let bestDistSq = Infinity;
            for (let rIdx = 0; rIdx < points.length; rIdx++) {
              const p = points[rIdx];
              const dSq = (vx - p.x) ** 2 + (vy - p.y) ** 2 + (vz - p.z) ** 2;
              if (dSq < bestDistSq) {
                bestDistSq = dSq;
                bestRes = segResidues[rIdx];
              }
            }
            const hex = getResidueColorHexNumber(bestRes, colorScheme, chainIdx);
            _color.setHex(hex);
            colors[v * 3] = _color.r;
            colors[v * 3 + 1] = _color.g;
            colors[v * 3 + 2] = _color.b;
          }

          tubeGeo.setAttribute('color', new THREE.BufferAttribute(colors, 3));

          const tubeMat = new THREE.MeshPhysicalMaterial({
            vertexColors: true,
            roughness: ssType === 'HELIX' ? 0.24 : ssType === 'SHEET' ? 0.36 : 0.58,
            metalness: ssType === 'HELIX' ? 0.18 : 0.08,
            clearcoat: ssType === 'HELIX' ? 0.55 : 0.2,
          });

          const tubeMesh = new THREE.Mesh(tubeGeo, tubeMat);
          tubeMesh.castShadow = !ctx.isSoftwareGl;
          tubeMesh.receiveShadow = !ctx.isSoftwareGl;
          ctx.proteinGroup.add(tubeMesh);

          // Add sculpted 3D Beveled Beta-Sheet Arrowhead at C-terminus of SHEET segments
          if (representation === 'cartoon-sidechain' && ssType === 'SHEET' && points.length >= 2) {
            const tipPt = points[points.length - 1];
            const prevPt = points[points.length - 2];
            _dir.subVectors(tipPt, prevPt).normalize();
            const arrowGeo = new THREE.ConeGeometry(1.05, 1.8, 4);
            arrowGeo.rotateX(Math.PI / 2);
            const lastRes = segResidues[segResidues.length - 1];
            const arrowHex = getResidueColorHexNumber(lastRes, colorScheme, chainIdx);
            const arrowMat = new THREE.MeshPhysicalMaterial({
              color: arrowHex,
              roughness: 0.3,
              metalness: 0.15,
              clearcoat: 0.35,
            });
            const arrowMesh = new THREE.Mesh(arrowGeo, arrowMat);
            arrowMesh.position.copy(tipPt);
            _quat.setFromUnitVectors(new THREE.Vector3(0, 0, 1), _dir);
            arrowMesh.quaternion.copy(_quat);
            ctx.proteinGroup.add(arrowMesh);
          }
        }

        segStart = segEnd + 1;
      }
    }

    // B. Instanced Heavy Atoms (THREE.InstancedMesh) — Raycastable on Layer 1 for real-time residue picking
    const resByKey = new Map<string, ResidueRecord>();
    for (const r of residues) resByKey.set(r.key, r);

    const visibleAtoms = atoms.filter((a) => {
      if (representation === 'backbone-trace') {
        return a.name === 'CA' || a.isHet || a.name === 'N' || a.name === 'O';
      }
      return true;
    });

    const sphereGeo = new THREE.SphereGeometry(1, 18, 14);
    const atomMat = new THREE.MeshPhysicalMaterial({
      color: 0xffffff,
      roughness: representation === 'spacefill' ? 0.32 : 0.25,
      metalness: 0.18,
      clearcoat: 0.4,
    });

    const atomInstancedMesh = new THREE.InstancedMesh(
      sphereGeo,
      atomMat,
      Math.max(1, visibleAtoms.length)
    );
    atomInstancedMesh.layers.enable(1); // Enable raycasting on layer 1
    atomInstancedMesh.castShadow = !ctx.isSoftwareGl;
    atomInstancedMesh.receiveShadow = !ctx.isSoftwareGl;

    for (let i = 0; i < visibleAtoms.length; i++) {
      const atom = visibleAtoms[i];
      getAtomWorldPos(atom, _vecA);

      let radius = 0.36;
      if (representation === 'spacefill') {
        radius = atom.element === 'C' ? 1.55 : atom.element === 'O' ? 1.42 : atom.element === 'N' ? 1.48 : 1.75;
      } else if (representation === 'ball-and-stick') {
        radius = atom.isHet ? 0.56 : atom.name === 'CA' ? 0.44 : 0.36;
      } else if (representation === 'cartoon-sidechain') {
        radius = atom.isHet ? 0.62 : atom.name === 'CA' ? 0.38 : 0.28;
      } else {
        radius = atom.isHet ? 0.55 : 0.32;
      }

      _scale.setScalar(radius);
      _quat.identity();
      _mat4.compose(_vecA, _quat, _scale);
      atomInstancedMesh.setMatrixAt(i, _mat4);

      // Determine atom color
      let hex = ELEMENT_COLORS[atom.element] ?? 0x94a3b8;
      if (colorScheme !== 'element' && (atom.name === 'CA' || representation === 'spacefill')) {
        const res = resByKey.get(`${atom.chainId}:${atom.resSeq}`);
        if (res) {
          const chainIdx = Math.max(0, model.chains.indexOf(atom.chainId));
          hex = getResidueColorHexNumber(res, colorScheme, chainIdx);
        }
      }
      _color.setHex(hex);
      atomInstancedMesh.setColorAt(i, _color);
    }

    atomInstancedMesh.instanceMatrix.needsUpdate = true;
    if (atomInstancedMesh.instanceColor) atomInstancedMesh.instanceColor.needsUpdate = true;
    atomInstancedMesh.computeBoundingSphere();

    ctx.proteinGroup.add(atomInstancedMesh);
    ctx.atomInstancedMesh = atomInstancedMesh;
    ctx.instanceToAtomMap = visibleAtoms;

    // C. Instanced Covalent & Disulfide Bonds (THREE.InstancedMesh)
    if (representation === 'ball-and-stick' || representation === 'cartoon-sidechain') {
      const cylGeo = new THREE.CylinderGeometry(1, 1, 1, 10, 1, false);
      const bondMat = new THREE.MeshStandardMaterial({
        color: 0xffffff,
        roughness: 0.45,
        metalness: 0.35,
      });

      const bondCount = model.covalentBonds.length;
      if (bondCount > 0) {
        const bondInstanced = new THREE.InstancedMesh(cylGeo, bondMat, bondCount);
        for (let i = 0; i < bondCount; i++) {
          const b = model.covalentBonds[i];
          const atomA = atoms[b.atomIndexA];
          const atomB = atoms[b.atomIndexB];
          getAtomWorldPos(atomA, _vecA);
          getAtomWorldPos(atomB, _vecB);

          _mid.addVectors(_vecA, _vecB).multiplyScalar(0.5);
          _dir.subVectors(_vecB, _vecA);
          const len = _dir.length();
          _dir.normalize();
          _quat.setFromUnitVectors(_up, _dir);

          const bondRadius =
            b.order === 'disulfide'
              ? 0.18
              : representation === 'ball-and-stick'
              ? 0.14
              : 0.095;
          _scale.set(bondRadius, len, bondRadius);
          _mat4.compose(_mid, _quat, _scale);
          bondInstanced.setMatrixAt(i, _mat4);

          _color.setHex(b.order === 'disulfide' ? 0xfacc15 : 0x64748b);
          bondInstanced.setColorAt(i, _color);
        }
        bondInstanced.instanceMatrix.needsUpdate = true;
        if (bondInstanced.instanceColor) bondInstanced.instanceColor.needsUpdate = true;
        bondInstanced.computeBoundingSphere();
        ctx.proteinGroup.add(bondInstanced);
      }
    }

    ctx.requestRender();
  }, [model, representation, colorScheme, explodedScale]);

  // 4. Rebuild Hydrogen Bond Network 3D Visualization when filteredHydrogenBonds, showHydrogenBonds, or activeHBondId changes
  useEffect(() => {
    const ctx = threeCtxRef.current;
    if (!ctx) return;

    disposeGroupChildren(ctx.hbondGroup);

    if (!showHydrogenBonds || filteredHydrogenBonds.length === 0) {
      ctx.requestRender();
      return;
    }

    const cx = ctx.centroid.x;
    const cy = ctx.centroid.y;
    const cz = ctx.centroid.z;

    // Render each H-bond as 5 beaded dashed segments + Donor/Acceptor orbital rings
    const beadsPerBond = 5;
    const totalBeads = filteredHydrogenBonds.length * beadsPerBond;
    const beadGeo = new THREE.CylinderGeometry(0.11, 0.11, 0.28, 10);
    const beadMat = new THREE.MeshStandardMaterial({
      color: 0xffffff,
      roughness: 0.2,
      metalness: 0.15,
      emissive: 0x0f172a,
    });

    const beadInstanced = new THREE.InstancedMesh(beadGeo, beadMat, Math.max(1, totalBeads));
    let beadIdx = 0;

    const ringGeo = new THREE.TorusGeometry(0.42, 0.06, 10, 24);

    for (const hb of filteredHydrogenBonds) {
      const isHighlighted =
        hb.id === activeHBondId ||
        hb.donorResKey === selectedResidueKey ||
        hb.acceptorResKey === selectedResidueKey;

      _vecA.set(hb.donorAtom.x - cx, hb.donorAtom.y - cy, hb.donorAtom.z - cz);
      _vecB.set(hb.acceptorAtom.x - cx, hb.acceptorAtom.y - cy, hb.acceptorAtom.z - cz);

      _dir.subVectors(_vecB, _vecA);
      const dist = _dir.length();
      _dir.normalize();
      _quat.setFromUnitVectors(_up, _dir);

      const catHex = isHighlighted
        ? 0xfef08a
        : HBOND_CATEGORY_COLORS[hb.category] ?? 0x38bdf8;
      _color.setHex(catHex);

      const radiusScale = isHighlighted ? 1.55 : 1.0;

      for (let s = 0; s < beadsPerBond; s++) {
        const t = (s + 0.5) / beadsPerBond;
        _mid.lerpVectors(_vecA, _vecB, t);
        _scale.set(radiusScale, (dist / beadsPerBond) * 1.65, radiusScale);
        _mat4.compose(_mid, _quat, _scale);
        beadInstanced.setMatrixAt(beadIdx, _mat4);
        beadInstanced.setColorAt(beadIdx, _color);
        beadIdx++;
      }

      // Add subtle donor/acceptor orbital halos for highlighted or active H-bonds
      if (isHighlighted) {
        const haloMat = new THREE.MeshBasicMaterial({
          color: catHex,
          transparent: true,
          opacity: 0.85,
        });
        const donorHalo = new THREE.Mesh(ringGeo.clone(), haloMat);
        donorHalo.position.copy(_vecA);
        donorHalo.quaternion.copy(_quat);
        ctx.hbondGroup.add(donorHalo);

        const acceptorHalo = new THREE.Mesh(ringGeo.clone(), haloMat.clone());
        acceptorHalo.position.copy(_vecB);
        acceptorHalo.quaternion.copy(_quat);
        ctx.hbondGroup.add(acceptorHalo);
      }
    }

    beadInstanced.instanceMatrix.needsUpdate = true;
    if (beadInstanced.instanceColor) beadInstanced.instanceColor.needsUpdate = true;
    beadInstanced.computeBoundingSphere();
    ctx.hbondGroup.add(beadInstanced);

    ctx.requestRender();
  }, [model, showHydrogenBonds, filteredHydrogenBonds, selectedResidueKey, activeHBondId]);

  // 5. Rebuild Selected / Hovered Residue 3D Analytical Beacon (`THREE.Group` of 5 distinct parts)
  useEffect(() => {
    const ctx = threeCtxRef.current;
    if (!ctx) return;

    disposeGroupChildren(ctx.beaconGroup);

    const targetKey = selectedResidueKey || hoveredResidueKey;
    if (!targetKey) {
      ctx.requestRender();
      return;
    }

    const res = model.residues.find((r) => r.key === targetKey);
    if (!res) {
      ctx.requestRender();
      return;
    }

    const cx = ctx.centroid.x;
    const cy = ctx.centroid.y;
    const cz = ctx.centroid.z;

    const caPos = new THREE.Vector3(
      res.caAtom.x - cx,
      res.caAtom.y - cy,
      res.caAtom.z - cz
    );

    const isSelected = targetKey === selectedResidueKey;
    const accentHex = isSelected ? 0x38bdf8 : 0xfbbf24;

    // Part 1: Outer Azimuth Gimbal Ring
    const outerRingGeo = new THREE.TorusGeometry(2.15, 0.065, 12, 48);
    const outerRingMat = new THREE.MeshStandardMaterial({
      color: accentHex,
      emissive: accentHex,
      emissiveIntensity: 0.45,
      roughness: 0.2,
      metalness: 0.8,
    });
    const outerRing = new THREE.Mesh(outerRingGeo, outerRingMat);
    outerRing.position.copy(caPos);
    ctx.beaconGroup.add(outerRing);

    // Part 2: Orthogonal Meridian Gimbal Ring
    const innerRingGeo = new THREE.TorusGeometry(1.75, 0.05, 12, 48);
    innerRingGeo.rotateX(Math.PI / 2);
    const innerRingMat = new THREE.MeshStandardMaterial({
      color: 0xf8fafc,
      roughness: 0.4,
      metalness: 0.5,
    });
    const innerRing = new THREE.Mesh(innerRingGeo, innerRingMat);
    innerRing.position.copy(caPos);
    ctx.beaconGroup.add(innerRing);

    // Part 3: Analytical Coordinate Spire Pointer above the residue
    const spireGeo = new THREE.ConeGeometry(0.28, 0.9, 6);
    spireGeo.rotateX(Math.PI);
    const spireMat = new THREE.MeshStandardMaterial({
      color: accentHex,
      roughness: 0.25,
      metalness: 0.75,
    });
    const spire = new THREE.Mesh(spireGeo, spireMat);
    spire.position.set(caPos.x, caPos.y + 2.75, caPos.z);
    ctx.beaconGroup.add(spire);

    // Part 4 & 5: Highlighted All-Atom Sidechain Ball-and-Stick Cluster for the Selected Residue
    const atomGeo = new THREE.SphereGeometry(0.46, 16, 12);
    for (const atom of res.atoms) {
      const aMat = new THREE.MeshPhysicalMaterial({
        color: ELEMENT_COLORS[atom.element] ?? 0x38bdf8,
        emissive: accentHex,
        emissiveIntensity: 0.22,
        roughness: 0.15,
        metalness: 0.2,
        clearcoat: 0.8,
      });
      const aMesh = new THREE.Mesh(atomGeo, aMat);
      aMesh.position.set(atom.x - cx, atom.y - cy, atom.z - cz);
      ctx.beaconGroup.add(aMesh);
    }

    ctx.requestRender();
  }, [model, selectedResidueKey, hoveredResidueKey]);

  // 6. Smooth Damped Camera Choreography on residue selection or camera preset click
  useEffect(() => {
    const ctx = threeCtxRef.current;
    if (!ctx) return;

    const r = Math.max(14, ctx.boundingRadius);
    const preset = cameraPresetTrigger.preset;

    if (preset === 'overview') {
      ctx.desiredCamTarget = new THREE.Vector3(0, 0, 0);
      ctx.desiredCamPos = new THREE.Vector3(0, r * 0.35, r * 2.15);
    } else if (preset === 'top') {
      ctx.desiredCamTarget = new THREE.Vector3(0, 0, 0);
      ctx.desiredCamPos = new THREE.Vector3(0, r * 2.35, 0.01);
    } else if (preset === 'side') {
      ctx.desiredCamTarget = new THREE.Vector3(0, 0, 0);
      ctx.desiredCamPos = new THREE.Vector3(r * 2.2, r * 0.2, 0);
    } else if (preset === 'focus-selected' && selectedResidueKey) {
      const res = model.residues.find((item) => item.key === selectedResidueKey);
      if (res) {
        const tx = res.caAtom.x - ctx.centroid.x;
        const ty = res.caAtom.y - ctx.centroid.y;
        const tz = res.caAtom.z - ctx.centroid.z;
        ctx.desiredCamTarget = new THREE.Vector3(tx, ty, tz);
        // Position camera radially outward from centroid through residue
        const out = new THREE.Vector3(tx, ty + 2.5, tz + 6.0);
        if (out.lengthSq() < 1) out.set(0, 3, 12);
        out.normalize().multiplyScalar(14.5);
        ctx.desiredCamPos = new THREE.Vector3(tx + out.x, ty + out.y, tz + out.z);
      }
    } else if (preset === 'focus-hbonds' && filteredHydrogenBonds.length > 0) {
      const hb =
        filteredHydrogenBonds.find((b) => b.id === activeHBondId) ||
        filteredHydrogenBonds[0];
      const tx = (hb.donorAtom.x + hb.acceptorAtom.x) * 0.5 - ctx.centroid.x;
      const ty = (hb.donorAtom.y + hb.acceptorAtom.y) * 0.5 - ctx.centroid.y;
      const tz = (hb.donorAtom.z + hb.acceptorAtom.z) * 0.5 - ctx.centroid.z;
      ctx.desiredCamTarget = new THREE.Vector3(tx, ty, tz);
      ctx.desiredCamPos = new THREE.Vector3(tx + 4.5, ty + 4.0, tz + 11.5);
    }

    ctx.requestRender();
  }, [cameraPresetTrigger, model, selectedResidueKey, activeHBondId, filteredHydrogenBonds]);

  return (
    <div
      ref={canvasMountRef}
      data-testid="protein-webgl-canvas"
      data-model-id={model.id}
      data-selected-residue={selectedResidueKey || 'none'}
      data-hbond-count={filteredHydrogenBonds.length}
      className="fixed inset-0 w-screen h-screen overflow-hidden z-0"
      aria-label={`Interactive 3D molecular viewport displaying ${model.title}`}
    />
  );
};

function StagePositionDefault(stageGroup: THREE.Group) {
  stageGroup.position.set(0, -11.5, 0);
}

function disposeGroupChildren(group: THREE.Group) {
  while (group.children.length > 0) {
    const child = group.children.pop()!;
    if ((child as THREE.Mesh).geometry) {
      (child as THREE.Mesh).geometry.dispose();
    }
    const mat = (child as THREE.Mesh).material;
    if (Array.isArray(mat)) {
      mat.forEach((m) => m.dispose());
    } else if (mat) {
      mat.dispose();
    }
  }
}
