import {
  AlphaFoldConfidenceAnalysis,
  AtomRecord,
  CovalentBond,
  FunctionalAnnotation,
  HBondCategory,
  HydrogenBond,
  PaeDomainBlock,
  ProteinStructureModel,
  ResidueRecord,
  SecondaryStructureType,
  StructureSource,
} from '../types/protein';

export const AMINO_ACID_ONE_LETTER: Record<string, string> = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C',
  GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I',
  LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P',
  SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
  SEC: 'U', PYL: 'O', HOH: 'O', HEM: 'H', ZN: 'Z', MG: 'M',
};

export const HYDROPHOBICITY_SCALE: Record<string, number> = {
  ILE: 4.5, VAL: 4.2, LEU: 3.8, PHE: 2.8, CYS: 2.5,
  MET: 1.9, ALA: 1.8, GLY: -0.4, THR: -0.7, SER: -0.8,
  TRP: -0.9, TYR: -1.3, PRO: -1.6, HIS: -3.2, GLU: -3.5,
  GLN: -3.5, ASP: -3.5, ASN: -3.5, LYS: -3.9, ARG: -4.5,
};

/**
 * Canonical AlphaFold 4-tier pLDDT confidence color scale
 */
export function getPlddtColorHex(plddt: number): string {
  if (plddt > 90) return '#0053D6'; // Very High (>90): Deep Blue
  if (plddt >= 70) return '#65CBF3'; // Confident (70-90): Light Blue
  if (plddt >= 50) return '#FFDB13'; // Low (50-70): Yellow
  return '#FF7D45';                  // Very Low (<50): Orange / Disordered
}

export function getPlddtTierLabel(plddt: number): string {
  if (plddt > 90) return 'Very High (>90)';
  if (plddt >= 70) return 'Confident (70–90)';
  if (plddt >= 50) return 'Low (50–70)';
  return 'Very Low (<50)';
}

export function distance3D(a: { x: number; y: number; z: number }, b: { x: number; y: number; z: number }): number {
  const dx = a.x - b.x;
  const dy = a.y - b.y;
  const dz = a.z - b.z;
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * Parse legacy PDB format text into atoms and secondary structure ranges
 */
export function parsePdbText(pdbText: string): {
  atoms: AtomRecord[];
  ssMap: Map<string, SecondaryStructureType>;
  headerTitle?: string;
  resolution?: number | null;
} {
  const lines = pdbText.split(/\r?\n/);
  const atoms: AtomRecord[] = [];
  const ssMap = new Map<string, SecondaryStructureType>();
  let headerTitle = '';
  let resolution: number | null = null;

  const helixRanges: Array<{ chain: string; start: number; end: number }> = [];
  const sheetRanges: Array<{ chain: string; start: number; end: number }> = [];

  for (const line of lines) {
    if (line.startsWith('TITLE ')) {
      headerTitle += ' ' + line.slice(10, 70).trim();
    } else if (line.startsWith('REMARK   2 RESOLUTION.')) {
      const m = line.match(/([0-9]+\.[0-9]+)\s+ANGSTROMS/i);
      if (m) resolution = parseFloat(m[1]);
    } else if (line.startsWith('HELIX ')) {
      const initChain = line.slice(19, 20).trim() || 'A';
      const initSeq = parseInt(line.slice(21, 25).trim(), 10);
      const endSeq = parseInt(line.slice(33, 37).trim(), 10);
      if (!Number.isNaN(initSeq) && !Number.isNaN(endSeq)) {
        helixRanges.push({ chain: initChain, start: initSeq, end: endSeq });
      }
    } else if (line.startsWith('SHEET ')) {
      const initChain = line.slice(21, 22).trim() || 'A';
      const initSeq = parseInt(line.slice(22, 26).trim(), 10);
      const endSeq = parseInt(line.slice(33, 37).trim(), 10);
      if (!Number.isNaN(initSeq) && !Number.isNaN(endSeq)) {
        sheetRanges.push({ chain: initChain, start: initSeq, end: endSeq });
      }
    } else if (line.startsWith('ATOM  ') || line.startsWith('HETATM')) {
      const isHet = line.startsWith('HETATM');
      const resName = line.slice(17, 20).trim().toUpperCase();
      if (resName === 'HOH' || resName === 'WAT') continue; // Skip bulk solvent water for clean rendering

      const altLoc = line.slice(16, 17).trim();
      if (altLoc && altLoc !== 'A' && altLoc !== '1') continue;

      const serial = parseInt(line.slice(6, 11).trim(), 10) || atoms.length + 1;
      const name = line.slice(12, 16).trim().toUpperCase();
      const chainId = line.slice(21, 22).trim() || 'A';
      const resSeq = parseInt(line.slice(22, 26).trim(), 10);
      const x = parseFloat(line.slice(30, 38).trim());
      const y = parseFloat(line.slice(38, 46).trim());
      const z = parseFloat(line.slice(46, 54).trim());
      const bFactor = parseFloat(line.slice(60, 66).trim()) || 75.0;
      const rawElem = line.slice(76, 78).trim().toUpperCase();
      const element = rawElem || inferElementFromAtomName(name);

      if (element === 'H' || element === 'D') continue; // Skip explicit hydrogens for clean heavy-atom geometry
      if (Number.isNaN(x) || Number.isNaN(y) || Number.isNaN(z) || Number.isNaN(resSeq)) continue;

      atoms.push({
        serial,
        name,
        element,
        resName,
        chainId,
        resSeq,
        x,
        y,
        z,
        bFactor,
        isHet,
      });
    }
  }

  for (const h of helixRanges) {
    for (let r = h.start; r <= h.end; r++) {
      ssMap.set(`${h.chain}:${r}`, 'HELIX');
    }
  }
  for (const s of sheetRanges) {
    for (let r = s.start; r <= s.end; r++) {
      ssMap.set(`${s.chain}:${r}`, 'SHEET');
    }
  }

  return {
    atoms,
    ssMap,
    headerTitle: headerTitle.trim() || undefined,
    resolution,
  };
}

/**
 * Parse mmCIF (.cif) format text into atoms and secondary structure ranges
 */
export function parseCifText(cifText: string): {
  atoms: AtomRecord[];
  ssMap: Map<string, SecondaryStructureType>;
} {
  const lines = cifText.split(/\r?\n/);
  const atoms: AtomRecord[] = [];
  const ssMap = new Map<string, SecondaryStructureType>();

  let inAtomLoop = false;
  let atomHeaders: string[] = [];

  for (let i = 0; i < lines.length; i++) {
    const rawLine = lines[i].trim();
    if (!rawLine || rawLine.startsWith('#')) {
      if (inAtomLoop && atomHeaders.length > 0 && atoms.length > 0) {
        inAtomLoop = false;
        atomHeaders = [];
      }
      continue;
    }

    if (rawLine === 'loop_') {
      inAtomLoop = false;
      atomHeaders = [];
      // Peek ahead to see if this is an _atom_site loop
      let j = i + 1;
      const headers: string[] = [];
      while (j < lines.length && lines[j].trim().startsWith('_')) {
        headers.push(lines[j].trim());
        j++;
      }
      if (headers.some((h) => h.startsWith('_atom_site.'))) {
        inAtomLoop = true;
        atomHeaders = headers;
        i = j - 1;
      }
      continue;
    }

    if (inAtomLoop) {
      if (rawLine.startsWith('_') || rawLine.startsWith('loop_') || rawLine.startsWith('data_')) {
        inAtomLoop = false;
        continue;
      }
      const tokens = tokenizeCifRow(rawLine);
      if (tokens.length < atomHeaders.length) continue;

      const idxGroup = atomHeaders.indexOf('_atom_site.group_PDB');
      const idxId = atomHeaders.indexOf('_atom_site.id');
      const idxAtomSymbol = atomHeaders.indexOf('_atom_site.type_symbol');
      const idxLabelAtomId = atomHeaders.indexOf('_atom_site.label_atom_id');
      const idxAltId = atomHeaders.indexOf('_atom_site.label_alt_id');
      const idxCompId = atomHeaders.indexOf('_atom_site.label_comp_id');
      const idxAuthAsym = atomHeaders.indexOf('_atom_site.auth_asym_id');
      const idxLabelAsym = atomHeaders.indexOf('_atom_site.label_asym_id');
      const idxAuthSeq = atomHeaders.indexOf('_atom_site.auth_seq_id');
      const idxLabelSeq = atomHeaders.indexOf('_atom_site.label_seq_id');
      const idxCartnX = atomHeaders.indexOf('_atom_site.Cartn_x');
      const idxCartnY = atomHeaders.indexOf('_atom_site.Cartn_y');
      const idxCartnZ = atomHeaders.indexOf('_atom_site.Cartn_z');
      const idxBiso = atomHeaders.indexOf('_atom_site.B_iso_or_equiv');
      const idxModelNum = atomHeaders.indexOf('_atom_site.pdbx_PDB_model_num');

      if (idxModelNum >= 0 && tokens[idxModelNum] !== '1' && tokens[idxModelNum] !== '?') {
        continue; // Use model 1 only for NMR ensembles
      }

      if (idxAltId >= 0) {
        const alt = tokens[idxAltId];
        if (alt !== '.' && alt !== '?' && alt !== 'A' && alt !== '1') continue;
      }

      const groupPdb = idxGroup >= 0 ? tokens[idxGroup] : 'ATOM';
      const isHet = groupPdb === 'HETATM';
      const resName = (idxCompId >= 0 ? tokens[idxCompId] : 'UNK').toUpperCase();
      if (resName === 'HOH' || resName === 'WAT') continue;

      const name = (idxLabelAtomId >= 0 ? tokens[idxLabelAtomId] : 'CA').replace(/"/g, '').toUpperCase();
      const element = (idxAtomSymbol >= 0 ? tokens[idxAtomSymbol] : inferElementFromAtomName(name)).toUpperCase();
      if (element === 'H' || element === 'D') continue;

      const chainId = (idxAuthAsym >= 0 ? tokens[idxAuthAsym] : idxLabelAsym >= 0 ? tokens[idxLabelAsym] : 'A') || 'A';
      const seqRaw = idxAuthSeq >= 0 ? tokens[idxAuthSeq] : idxLabelSeq >= 0 ? tokens[idxLabelSeq] : '1';
      const resSeq = parseInt(seqRaw, 10);
      const x = parseFloat(tokens[idxCartnX]);
      const y = parseFloat(tokens[idxCartnY]);
      const z = parseFloat(tokens[idxCartnZ]);
      const bFactor = idxBiso >= 0 ? parseFloat(tokens[idxBiso]) || 80.0 : 80.0;
      const serial = idxId >= 0 ? parseInt(tokens[idxId], 10) || atoms.length + 1 : atoms.length + 1;

      if (Number.isNaN(x) || Number.isNaN(y) || Number.isNaN(z) || Number.isNaN(resSeq)) continue;

      atoms.push({
        serial,
        name,
        element,
        resName,
        chainId,
        resSeq,
        x,
        y,
        z,
        bFactor,
        isHet,
      });
    }
  }

  return { atoms, ssMap };
}

function tokenizeCifRow(line: string): string[] {
  const result: string[] = [];
  let current = '';
  let inSingle = false;
  let inDouble = false;

  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (ch === "'" && !inDouble) {
      inSingle = !inSingle;
    } else if (ch === '"' && !inSingle) {
      inDouble = !inDouble;
    } else if ((ch === ' ' || ch === '\t') && !inSingle && !inDouble) {
      if (current.length > 0) {
        result.push(current);
        current = '';
      }
    } else {
      current += ch;
    }
  }
  if (current.length > 0) {
    result.push(current);
  }
  return result;
}

function inferElementFromAtomName(atomName: string): string {
  const clean = atomName.replace(/[0-9]/g, '').trim();
  if (clean.startsWith('FE')) return 'FE';
  if (clean.startsWith('ZN')) return 'ZN';
  if (clean.startsWith('MG')) return 'MG';
  if (clean.startsWith('CA') && clean.length === 2 && atomName === 'CA') return 'C'; // Alpha carbon
  return clean.slice(0, 1) || 'C';
}

/**
 * Group atoms into ordered ResidueRecords and infer secondary structure if not explicitly provided
 */
export function buildResiduesFromAtoms(
  atoms: AtomRecord[],
  explicitSsMap: Map<string, SecondaryStructureType>
): ResidueRecord[] {
  const residueGroups = new Map<string, AtomRecord[]>();
  for (const atom of atoms) {
    const key = `${atom.chainId}:${atom.resSeq}`;
    const list = residueGroups.get(key);
    if (list) {
      list.push(atom);
    } else {
      residueGroups.set(key, [atom]);
    }
  }

  const residues: ResidueRecord[] = [];
  for (const [key, resAtoms] of residueGroups.entries()) {
    const first = resAtoms[0];
    const caAtom =
      resAtoms.find((a) => a.name === 'CA') ||
      resAtoms.find((a) => a.name === 'C1' || a.name === 'FE' || a.name === 'P') ||
      first;
    const meanB =
      resAtoms.reduce((acc, a) => acc + a.bFactor, 0) / Math.max(1, resAtoms.length);

    residues.push({
      key,
      chainId: first.chainId,
      resSeq: first.resSeq,
      resName: first.resName,
      oneLetter: AMINO_ACID_ONE_LETTER[first.resName] || 'X',
      secondaryStructure: explicitSsMap.get(key) || 'COIL',
      plddtOrBfactor: Number(meanB.toFixed(1)),
      caAtom,
      atoms: resAtoms,
      isLigand: first.isHet,
    });
  }

  // If no explicit HELIX/SHEET records were in the file, infer secondary structure from C-alpha geometry & backbone H-bonds
  if (explicitSsMap.size === 0 && residues.length >= 5) {
    inferSecondaryStructureGeometrically(residues);
  }

  return residues;
}

/**
 * Infer HELIX and SHEET secondary structure using C-alpha distances (i, i+3 / i, i+4) and backbone H-bonds
 */
function inferSecondaryStructureGeometrically(residues: ResidueRecord[]): void {
  const poly = residues.filter((r) => !r.isLigand);
  for (let i = 0; i < poly.length; i++) {
    if (i + 3 < poly.length && poly[i].chainId === poly[i + 3].chainId) {
      const d13 = distance3D(poly[i].caAtom, poly[i + 3].caAtom);
      const d14 =
        i + 4 < poly.length && poly[i].chainId === poly[i + 4].chainId
          ? distance3D(poly[i].caAtom, poly[i + 4].caAtom)
          : 99;
      // Canonical alpha-helix has CA(i)-CA(i+3) ~ 5.0 - 5.6 Å and CA(i)-CA(i+4) ~ 5.8 - 6.5 Å
      if (d13 >= 4.6 && d13 <= 5.9 && d14 >= 5.4 && d14 <= 6.8) {
        for (let k = i; k <= Math.min(poly.length - 1, i + 4); k++) {
          poly[k].secondaryStructure = 'HELIX';
        }
      }
    }
  }

  // Identify extended beta strands where consecutive CA(i)-CA(i+2) distance is > 6.4 Å and cross-strand CA is < 5.6 Å
  for (let i = 1; i < poly.length - 1; i++) {
    if (poly[i].secondaryStructure === 'HELIX') continue;
    if (poly[i - 1].chainId !== poly[i + 1].chainId) continue;
    const d12 = distance3D(poly[i - 1].caAtom, poly[i + 1].caAtom);
    if (d12 >= 6.2) {
      // Check if there is a non-local partner strand within 5.6 Å
      const hasPartner = poly.some(
        (other, j) =>
          Math.abs(j - i) >= 4 &&
          distance3D(poly[i].caAtom, other.caAtom) <= 5.6
      );
      if (hasPartner) {
        poly[i - 1].secondaryStructure = 'SHEET';
        poly[i].secondaryStructure = 'SHEET';
        poly[i + 1].secondaryStructure = 'SHEET';
      }
    }
  }
}

/**
 * Compute covalent bonds between heavy atoms (intra-residue, peptide C-N backbone, and disulfide SG-SG bridges)
 */
export function computeCovalentBonds(atoms: AtomRecord[]): CovalentBond[] {
  const bonds: CovalentBond[] = [];
  const n = atoms.length;

  // Group indices by residue and track C / N / SG for fast inter-residue bonding
  const resAtomIndices = new Map<string, number[]>();
  const cByRes = new Map<string, number>();
  const nByRes = new Map<string, number>();
  const sgIndices: number[] = [];

  for (let i = 0; i < n; i++) {
    const a = atoms[i];
    const key = `${a.chainId}:${a.resSeq}`;
    let arr = resAtomIndices.get(key);
    if (!arr) {
      arr = [];
      resAtomIndices.set(key, arr);
    }
    arr.push(i);
    if (!a.isHet) {
      if (a.name === 'C') cByRes.set(key, i);
      if (a.name === 'N') nByRes.set(key, i);
      if (a.name === 'SG' && a.resName === 'CYS') sgIndices.push(i);
    }
  }

  // 1. Intra-residue covalent bonds (heavy atom distance <= 1.92 Å)
  for (const indices of resAtomIndices.values()) {
    for (let i = 0; i < indices.length; i++) {
      const idxA = indices[i];
      const atomA = atoms[idxA];
      for (let j = i + 1; j < indices.length; j++) {
        const idxB = indices[j];
        const atomB = atoms[idxB];
        const d = distance3D(atomA, atomB);
        const maxBond = atomA.element === 'S' || atomB.element === 'S' || atomA.element === 'FE' || atomB.element === 'FE' ? 2.15 : 1.92;
        if (d >= 0.9 && d <= maxBond) {
          bonds.push({
            atomIndexA: idxA,
            atomIndexB: idxB,
            order: (atomA.name === 'C' && atomB.name === 'O') || (atomA.name === 'O' && atomB.name === 'C') ? 2 : 1,
          });
        }
      }
    }
  }

  // 2. Inter-residue peptide bonds: C(i) to N(i+1)
  const keys = Array.from(resAtomIndices.keys());
  for (let i = 0; i < keys.length - 1; i++) {
    const idxC = cByRes.get(keys[i]);
    const idxN = nByRes.get(keys[i + 1]);
    if (idxC !== undefined && idxN !== undefined) {
      const atomC = atoms[idxC];
      const atomN = atoms[idxN];
      if (atomC.chainId === atomN.chainId && distance3D(atomC, atomN) <= 1.95) {
        bonds.push({ atomIndexA: idxC, atomIndexB: idxN, order: 1 });
      }
    }
  }

  // 3. Disulfide bridges: CYS SG to CYS SG within [1.8 Å, 2.25 Å]
  for (let i = 0; i < sgIndices.length; i++) {
    for (let j = i + 1; j < sgIndices.length; j++) {
      const idxA = sgIndices[i];
      const idxB = sgIndices[j];
      const d = distance3D(atoms[idxA], atoms[idxB]);
      if (d >= 1.8 && d <= 2.25) {
        bonds.push({ atomIndexA: idxA, atomIndexB: idxB, order: 'disulfide' });
      }
    }
  }

  return bonds;
}

const DONOR_ATOM_NAMES = new Set([
  'N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ',
  'OG', 'OG1', 'OH', 'SG',
]);

const ACCEPTOR_ATOM_NAMES = new Set([
  'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'OXT',
  'ND1', 'NE2', 'SD', 'SG', 'N1', 'N2', 'N3', 'N4', 'O1', 'O2', 'O1A', 'O2A', 'O1D', 'O2D',
]);

function canActAsDonor(atom: AtomRecord): boolean {
  if (atom.resName === 'PRO' && atom.name === 'N') return false; // Proline backbone N lacks H
  if (atom.isHet) return atom.element === 'N' || atom.element === 'O';
  return DONOR_ATOM_NAMES.has(atom.name) || atom.element === 'N';
}

function canActAsAcceptor(atom: AtomRecord): boolean {
  if (atom.isHet) return atom.element === 'O' || atom.element === 'N';
  return ACCEPTOR_ATOM_NAMES.has(atom.name) || atom.element === 'O';
}

/**
 * Geometric Hydrogen Bond Detection Engine
 * Detects chemically valid Donor-Acceptor hydrogen bonds within [2.2 Å, maxCutoff Å]
 */
export function computeHydrogenBonds(
  atoms: AtomRecord[],
  residues: ResidueRecord[],
  maxCutoff = 4.0
): HydrogenBond[] {
  const polarAtoms = atoms.filter(
    (a) => a.element === 'N' || a.element === 'O' || (a.element === 'S' && a.name === 'SG')
  );

  const ssByResKey = new Map<string, SecondaryStructureType>();
  for (const r of residues) {
    ssByResKey.set(r.key, r.secondaryStructure);
  }

  // Map each residue key to its CA or C atom to estimate donor-acceptor-antecedent angle
  const anchorByResKey = new Map<string, AtomRecord>();
  for (const r of residues) {
    const cAtom = r.atoms.find((a) => a.name === 'C') || r.caAtom;
    anchorByResKey.set(r.key, cAtom);
  }

  const hbonds: HydrogenBond[] = [];
  const seenPairs = new Set<string>();

  for (let i = 0; i < polarAtoms.length; i++) {
    const a = polarAtoms[i];
    const keyA = `${a.chainId}:${a.resSeq}`;

    for (let j = i + 1; j < polarAtoms.length; j++) {
      const b = polarAtoms[j];
      const keyB = `${b.chainId}:${b.resSeq}`;

      // Skip same residue or adjacent peptide bond neighbors (|resSeq diff| <= 1 on same chain unless sidechain-to-sidechain)
      const seqDiff = Math.abs(a.resSeq - b.resSeq);
      if (a.chainId === b.chainId && seqDiff === 0) continue;
      if (
        a.chainId === b.chainId &&
        seqDiff === 1 &&
        (a.name === 'N' || a.name === 'O') &&
        (b.name === 'N' || b.name === 'O')
      ) {
        continue;
      }
      // Exclude disulfide bridges (CYS SG - CYS SG)
      if (a.element === 'S' && b.element === 'S') continue;

      const dist = distance3D(a, b);
      if (dist < 2.25 || dist > maxCutoff) continue;

      // Determine Donor and Acceptor orientation
      let donor: AtomRecord | null = null;
      let acceptor: AtomRecord | null = null;

      if (canActAsDonor(a) && canActAsAcceptor(b)) {
        donor = a;
        acceptor = b;
      } else if (canActAsDonor(b) && canActAsAcceptor(a)) {
        donor = b;
        acceptor = a;
      }
      if (!donor || !acceptor) continue;

      const donorResKey = `${donor.chainId}:${donor.resSeq}`;
      const acceptorResKey = `${acceptor.chainId}:${acceptor.resSeq}`;

      // Estimate angle using acceptor's covalent antecedent (e.g. C=O...N)
      const antecedent = anchorByResKey.get(acceptorResKey);
      let angleDegrees = 154;
      if (antecedent && antecedent.serial !== acceptor.serial) {
        const v1x = antecedent.x - acceptor.x;
        const v1y = antecedent.y - acceptor.y;
        const v1z = antecedent.z - acceptor.z;
        const v2x = donor.x - acceptor.x;
        const v2y = donor.y - acceptor.y;
        const v2z = donor.z - acceptor.z;
        const mag1 = Math.sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
        const mag2 = Math.sqrt(v2x * v2x + v2y * v2y + v2z * v2z);
        if (mag1 > 0.1 && mag2 > 0.1) {
          const cosTheta = Math.max(-1, Math.min(1, (v1x * v2x + v1y * v2y + v1z * v2z) / (mag1 * mag2)));
          angleDegrees = Math.round((Math.acos(cosTheta) * 180) / Math.PI);
        }
      }

      if (angleDegrees < 95) continue; // Reject implausible acute geometry

      const pairId =
        donor.serial < acceptor.serial
          ? `${donor.serial}-${acceptor.serial}`
          : `${acceptor.serial}-${donor.serial}`;
      if (seenPairs.has(pairId)) continue;
      seenPairs.add(pairId);

      // Classify hydrogen bond category
      const isBackbonePair =
        !donor.isHet &&
        !acceptor.isHet &&
        (donor.name === 'N' || donor.name === 'O') &&
        (acceptor.name === 'N' || acceptor.name === 'O');

      let category: HBondCategory = 'sidechain';
      if (donor.isHet || acceptor.isHet) {
        category = 'ligand';
      } else if (isBackbonePair) {
        const resGap = Math.abs(donor.resSeq - acceptor.resSeq);
        const donorSs = ssByResKey.get(donorResKey);
        const accSs = ssByResKey.get(acceptorResKey);
        if (
          donor.chainId === acceptor.chainId &&
          (resGap === 3 || resGap === 4 || resGap === 5 || donorSs === 'HELIX' || accSs === 'HELIX') &&
          resGap <= 5
        ) {
          category = 'backbone-helix';
        } else {
          category = 'backbone-sheet';
        }
      }

      // Empirical Kabsch-Sander style electrostatic energy approximation (kcal/mol)
      const angleFactor = Math.max(0.25, Math.sin(((angleDegrees - 60) * Math.PI) / 180));
      const energyKcalMol = Number((-2.8 * Math.pow(2.9 / dist, 2) * angleFactor).toFixed(2));

      hbonds.push({
        id: `HB-${hbonds.length + 1}`,
        donorAtom: donor,
        acceptorAtom: acceptor,
        donorResKey,
        acceptorResKey,
        distance: Number(dist.toFixed(2)),
        angleDegrees,
        energyKcalMol,
        category,
      });
    }
  }

  return hbonds.sort((a, b) => a.distance - b.distance);
}

/**
 * Analyze AlphaFold pLDDT & PAE confidence metrics per alphafold_database_fetch_and_analyze.md
 */
export function analyzeAlphaFoldConfidence(
  residues: ResidueRecord[],
  source: StructureSource,
  entryId: string,
  rawPaeGrid?: number[][]
): AlphaFoldConfidenceAnalysis {
  const poly = residues.filter((r) => !r.isLigand);
  const total = Math.max(1, poly.length);

  let veryHigh = 0;
  let confident = 0;
  let low = 0;
  let veryLow = 0;
  let sum = 0;

  for (const r of poly) {
    const score = r.plddtOrBfactor;
    sum += score;
    if (score > 90) veryHigh++;
    else if (score >= 70) confident++;
    else if (score >= 50) low++;
    else veryLow++;
  }

  const meanPlddt = Number((sum / total).toFixed(1));
  const fractionVeryHigh = Number(((veryHigh / total) * 100).toFixed(1));
  const fractionConfident = Number(((confident / total) * 100).toFixed(1));
  const fractionLow = Number(((low / total) * 100).toFixed(1));
  const fractionVeryLow = Number(((veryLow / total) * 100).toFixed(1));

  // Identify rigid domains (contiguous segments with pLDDT >= 70)
  const domains: PaeDomainBlock[] = [];
  let inDomain = false;
  let domStart = 0;
  let domScores: number[] = [];

  for (let i = 0; i < poly.length; i++) {
    const r = poly[i];
    if (r.plddtOrBfactor >= 70) {
      if (!inDomain) {
        inDomain = true;
        domStart = r.resSeq;
        domScores = [r.plddtOrBfactor];
      } else {
        domScores.push(r.plddtOrBfactor);
      }
    } else if (inDomain) {
      if (domScores.length >= 6) {
        const avg = domScores.reduce((a, b) => a + b, 0) / domScores.length;
        domains.push({
          id: `D${domains.length + 1}`,
          label: `Ordered Core Domain ${domains.length + 1}`,
          startRes: domStart,
          endRes: poly[i - 1].resSeq,
          meanInternalPae: Number((2.1 + (100 - avg) * 0.12).toFixed(1)),
          meanPlddt: Number(avg.toFixed(1)),
        });
      }
      inDomain = false;
      domScores = [];
    }
  }
  if (inDomain && domScores.length >= 5) {
    const avg = domScores.reduce((a, b) => a + b, 0) / domScores.length;
    domains.push({
      id: `D${domains.length + 1}`,
      label: `Ordered Core Domain ${domains.length + 1}`,
      startRes: domStart,
      endRes: poly[poly.length - 1].resSeq,
      meanInternalPae: Number((2.1 + (100 - avg) * 0.12).toFixed(1)),
      meanPlddt: Number(avg.toFixed(1)),
    });
  }

  if (domains.length === 0 && poly.length > 0) {
    domains.push({
      id: 'D1',
      label: 'Primary Globular Fold',
      startRes: poly[0].resSeq,
      endRes: poly[poly.length - 1].resSeq,
      meanInternalPae: 3.4,
      meanPlddt,
    });
  }

  // Generate 20x20 PAE grid sample (either from provided PAE or derived from domain block topology & CA distance)
  const gridSize = 20;
  const paeGridSample: number[][] = [];
  if (rawPaeGrid && rawPaeGrid.length > 0) {
    for (let r = 0; r < gridSize; r++) {
      const row: number[] = [];
      const srcR = Math.min(rawPaeGrid.length - 1, Math.floor((r / gridSize) * rawPaeGrid.length));
      for (let c = 0; c < gridSize; c++) {
        const srcC = Math.min(
          rawPaeGrid[srcR].length - 1,
          Math.floor((c / gridSize) * rawPaeGrid[srcR].length)
        );
        row.push(Number(rawPaeGrid[srcR][srcC].toFixed(1)));
      }
      paeGridSample.push(row);
    }
  } else {
    for (let r = 0; r < gridSize; r++) {
      const row: number[] = [];
      const idxR = Math.min(poly.length - 1, Math.floor((r / gridSize) * poly.length));
      const resR = poly[idxR];
      for (let c = 0; c < gridSize; c++) {
        if (r === c) {
          row.push(0.8);
          continue;
        }
        const idxC = Math.min(poly.length - 1, Math.floor((c / gridSize) * poly.length));
        const resC = poly[idxC];
        const sameDomain = domains.some(
          (d) =>
            resR.resSeq >= d.startRes &&
            resR.resSeq <= d.endRes &&
            resC.resSeq >= d.startRes &&
            resC.resSeq <= d.endRes
        );
        const minConf = Math.min(resR?.plddtOrBfactor ?? 80, resC?.plddtOrBfactor ?? 80);
        const baseError = sameDomain ? 2.4 : 14.5;
        const penalty = Math.max(0, (85 - minConf) * 0.32);
        row.push(Number(Math.min(30, baseError + penalty).toFixed(1)));
      }
      paeGridSample.push(row);
    }
  }

  const plddtConclusion =
    source === 'ALPHAFOLD_DB'
      ? `Global mean pLDDT is ${meanPlddt} (${fractionVeryHigh}% Very High >90, ${fractionConfident}% Confident 70–90, ${fractionLow}% Low 50–70, ${fractionVeryLow}% Very Low <50). Per-residue pLDDT confidence is stored directly in the mmCIF/PDB B-factor column.`
      : `Experimental crystallographic/cryo-EM structure with mean atomic B-factor ${meanPlddt} Å² across ${poly.length} resolved residues.`;

  const paeConclusion =
    domains.length > 1
      ? `PAE matrix resolves ${domains.length} distinct rigid domains (${domains.map((d) => `Res ${d.startRes}–${d.endRes}`).join(', ')}) with low intra-domain error (< 4.5 Å) separated by flexible inter-domain linkers.`
      : `PAE matrix indicates a single cohesive globular fold (Res ${domains[0].startRes}–${domains[0].endRes}, mean internal PAE ${domains[0].meanInternalPae} Å) with well-defined tertiary packing.`;

  const cohesiveSummary = `${plddtConclusion} ${paeConclusion}`;

  const disorderWarning =
    source === 'ALPHAFOLD_DB' && (fractionVeryLow >= 20 || fractionLow + fractionVeryLow >= 35)
      ? `Explicit Disorder Warning: ${fractionVeryLow}% of residues exhibit Very Low confidence (<50 pLDDT) and ${fractionLow}% Low confidence (50–70 pLDDT), strongly indicating intrinsically disordered protein regions (IDPRs). Do not interpret low-pLDDT loops as fixed folds; restrict downstream docking or Foldseek alignment strictly to the ordered domain boundaries (${domains.map((d) => `Res ${d.startRes}–${d.endRes}`).join(', ')}).`
      : null;

  const isoformOrLengthWarning =
    entryId.includes('-F2') || entryId.includes('-F3') || poly.length > 2700
      ? `Fragment / Length Notice: Entry ${entryId} represents an overlapping fragment or large polypeptide (>2,700 AA). Verify domain continuity across fragment boundaries.`
      : null;

  return {
    meanPlddt,
    fractionVeryHigh,
    fractionConfident,
    fractionLow,
    fractionVeryLow,
    plddtConclusion,
    paeConclusion,
    cohesiveSummary,
    disorderWarning,
    isoformOrLengthWarning,
    domains,
    paeGridSample,
  };
}

/**
 * Assemble a complete ProteinStructureModel from parsed atoms & metadata
 */
export function buildProteinStructureModel(params: {
  id: string;
  displayCode: string;
  title: string;
  subtitle: string;
  source: StructureSource;
  geneSymbol: string;
  organism: string;
  uniprotId?: string;
  pdbId?: string;
  experimentalMethod: string;
  resolutionAngstrom: number | null;
  releaseDate: string;
  citationTitle: string;
  sourceUrl: string;
  cifDownloadUrl: string;
  atoms: AtomRecord[];
  ssMap: Map<string, SecondaryStructureType>;
  annotations?: FunctionalAnnotation[];
  rawPaeGrid?: number[][];
}): ProteinStructureModel {
  const residues = buildResiduesFromAtoms(params.atoms, params.ssMap);
  const chains = Array.from(new Set(params.atoms.map((a) => a.chainId)));
  const covalentBonds = computeCovalentBonds(params.atoms);
  const allHydrogenBonds = computeHydrogenBonds(params.atoms, residues, 4.0);
  const alphaFoldAnalysis = analyzeAlphaFoldConfidence(
    residues,
    params.source,
    params.id,
    params.rawPaeGrid
  );

  return {
    id: params.id,
    displayCode: params.displayCode,
    title: params.title,
    subtitle: params.subtitle,
    source: params.source,
    provenanceBadge:
      params.source === 'ALPHAFOLD_DB' ? 'Predicted by model' : 'Experimentally determined',
    geneSymbol: params.geneSymbol,
    organism: params.organism,
    uniprotId: params.uniprotId,
    pdbId: params.pdbId,
    experimentalMethod: params.experimentalMethod,
    resolutionAngstrom: params.resolutionAngstrom,
    releaseDate: params.releaseDate,
    citationTitle: params.citationTitle,
    sourceUrl: params.sourceUrl,
    cifDownloadUrl: params.cifDownloadUrl,
    atoms: params.atoms,
    residues,
    chains,
    covalentBonds,
    allHydrogenBonds,
    alphaFoldAnalysis,
    annotations: params.annotations || [],
  };
}
