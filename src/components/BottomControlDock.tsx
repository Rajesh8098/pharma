import React, { useEffect, useRef, useState } from 'react';
import {
  ColorSchemeMode,
  HBondFilterMode,
  HydrogenBond,
  PdbSearchResultItem,
  ProteinStructureModel,
  RepresentationMode,
  ResidueRecord,
} from '../types/protein';
import { getPlddtColorHex } from '../utils/proteinParser';
import { CameraPresetType } from './ProteinViewport3D';
import { DockTabId } from './TopBar';

interface BottomControlDockProps {
  model: ProteinStructureModel;
  activeTab: DockTabId;
  dockCollapsed: boolean;
  representation: RepresentationMode;
  colorScheme: ColorSchemeMode;
  showHydrogenBonds: boolean;
  hbondCutoff: number;
  hbondFilter: HBondFilterMode;
  filteredHydrogenBonds: HydrogenBond[];
  selectedResidueKey: string | null;
  activeHBondId: string | null;
  lightBackground: boolean;
  explodedScale: number;
  queryExplanation: string;
  searchQuery: string;
  maxResolutionFilter: number;
  searchResults: PdbSearchResultItem[];
  searchExplanation: string;
  searchError: string | null;
  isSearching: boolean;
  isLoadingStructure: boolean;
  loadError: string | null;
  onSelectTab: (tab: DockTabId) => void;
  onToggleCollapse: () => void;
  onChangeRepresentation: (mode: RepresentationMode) => void;
  onChangeColorScheme: (scheme: ColorSchemeMode) => void;
  onToggleHydrogenBonds: () => void;
  onChangeHBondCutoff: (val: number) => void;
  onChangeHBondFilter: (filter: HBondFilterMode) => void;
  onSelectResidue: (residueKey: string) => void;
  onSelectHBond: (hbond: HydrogenBond) => void;
  onToggleLightBackground: () => void;
  onChangeExplodedScale: (val: number) => void;
  onTriggerCameraPreset: (preset: CameraPresetType) => void;
  onSearchQueryChange: (val: string) => void;
  onMaxResolutionChange: (val: number) => void;
  onSubmitSearch: (e: React.FormEvent) => void;
  onLoadAccession: (code: string, forceLive?: boolean) => void;
}

export const BottomControlDock: React.FC<BottomControlDockProps> = ({
  model,
  activeTab,
  dockCollapsed,
  representation,
  colorScheme,
  showHydrogenBonds,
  hbondCutoff,
  hbondFilter,
  filteredHydrogenBonds,
  selectedResidueKey,
  activeHBondId,
  lightBackground,
  explodedScale,
  queryExplanation,
  searchQuery,
  maxResolutionFilter,
  searchResults,
  searchExplanation,
  searchError,
  isSearching,
  isLoadingStructure,
  loadError,
  onSelectTab,
  onToggleCollapse,
  onChangeRepresentation,
  onChangeColorScheme,
  onToggleHydrogenBonds,
  onChangeHBondCutoff,
  onChangeHBondFilter,
  onSelectResidue,
  onSelectHBond,
  onToggleLightBackground,
  onChangeExplodedScale,
  onTriggerCameraPreset,
  onSearchQueryChange,
  onMaxResolutionChange,
  onSubmitSearch,
  onLoadAccession,
}) => {
  const [seqFilter, setSeqFilter] = useState('');
  const paeCanvasRef = useRef<HTMLCanvasElement | null>(null);

  const selectedResidue: ResidueRecord | undefined =
    model.residues.find((r) => r.key === selectedResidueKey) || model.residues[0];

  // Draw 20x20 PAE (Predicted Aligned Error) matrix onto HTML5 Canvas when confidence tab is open
  useEffect(() => {
    if (activeTab !== 'confidence' || dockCollapsed) return;
    const canvas = paeCanvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const grid = model.alphaFoldAnalysis.paeGridSample;
    const n = grid.length;
    const size = canvas.width;
    const cell = size / Math.max(1, n);

    ctx.clearRect(0, 0, size, size);
    for (let r = 0; r < n; r++) {
      for (let c = 0; c < n; c++) {
        const val = grid[r][c]; // 0 to 30 Å
        const t = Math.max(0, Math.min(1, val / 28));
        // Deep emerald/teal for low PAE (<5 Å) to light slate/amber for high PAE (>20 Å)
        const red = Math.round(13 + t * 215);
        const green = Math.round(148 - t * 90);
        const blue = Math.round(136 - t * 95);
        ctx.fillStyle = `rgb(${red}, ${green}, ${blue})`;
        ctx.fillRect(c * cell, r * cell, Math.ceil(cell), Math.ceil(cell));
      }
    }
  }, [activeTab, dockCollapsed, model]);

  const filteredResidues = model.residues.filter((r) => {
    if (!seqFilter.trim()) return true;
    const q = seqFilter.trim().toUpperCase();
    return (
      r.resName.includes(q) ||
      String(r.resSeq).includes(q) ||
      r.oneLetter === q ||
      r.secondaryStructure.includes(q)
    );
  });

  return (
    <section
      aria-label="Interactive structure and analysis dock"
      className="fixed bottom-2 sm:bottom-3 left-1/2 -translate-x-1/2 z-20 w-[min(1040px,95vw)] pointer-events-none flex flex-col items-center"
    >
      {/* Top Collapse / Expand Handle (Zen / Cinema Mode Drawer Toggle) */}
      <div className="pointer-events-auto flex items-center gap-1.5 mb-1">
        <button
          type="button"
          onClick={onToggleCollapse}
          className="px-3.5 py-1 text-xs font-medium text-slate-200 bg-[#0f172a]/95 hover:bg-slate-800 border border-slate-700 rounded-t-md transition-colors whitespace-nowrap shrink-0"
        >
          {dockCollapsed ? 'Show Controls (H)' : 'Hide Controls (H)'}
        </button>
      </div>

      {/* Main Unified Bottom Control Dock */}
      {!dockCollapsed && (
        <div className="pointer-events-auto w-full max-h-[38vh] sm:max-h-none overflow-y-auto bg-[#0f172a]/95 border border-slate-800 rounded-xl shadow-2xl">
          {/* Segmented Dock Tab Switcher + Camera Quick-Actions */}
          <div className="flex flex-wrap items-center justify-between gap-2 px-3 sm:px-4 py-2 bg-[#0b0f19]/90 border-b border-slate-800 sticky top-0 z-10">
            <div
              role="tablist"
              aria-label="Analysis modes"
              className="flex items-center gap-1 overflow-x-auto max-w-full"
            >
              {(
                [
                  { id: 'viewport', label: '3D Controls' },
                  {
                    id: 'sequence',
                    label: `Sequence (${model.residues.length})`,
                  },
                  {
                    id: 'hbonds',
                    label: `H-Bonds (${filteredHydrogenBonds.length})`,
                  },
                  { id: 'confidence', label: 'pLDDT & PAE' },
                  { id: 'search', label: 'PDB & AlphaFold Search' },
                ] as const
              ).map((tab) => {
                const isActive = activeTab === tab.id;
                return (
                  <button
                    key={tab.id}
                    role="tab"
                    aria-selected={isActive}
                    type="button"
                    onClick={() => onSelectTab(tab.id)}
                    className={`px-2.5 sm:px-3 py-1 sm:py-1.5 text-xs font-medium rounded-md transition-colors whitespace-nowrap shrink-0 ${
                      isActive
                        ? 'bg-cyan-500 text-slate-950 font-semibold'
                        : 'text-slate-300 hover:text-white hover:bg-slate-800'
                    }`}
                  >
                    {tab.label}
                  </button>
                );
              })}
            </div>

            {/* Camera & Stage Quick Bar */}
            <div className="hidden sm:flex items-center gap-1.5">
              <button
                type="button"
                onClick={() => onTriggerCameraPreset('overview')}
                className="px-2.5 py-1 text-[11px] font-medium text-slate-300 bg-slate-800 hover:bg-slate-700 rounded transition-colors whitespace-nowrap shrink-0"
              >
                Reset Camera
              </button>
              <button
                type="button"
                onClick={() => onTriggerCameraPreset('top')}
                className="px-2.5 py-1 text-[11px] font-medium text-slate-300 bg-slate-800 hover:bg-slate-700 rounded transition-colors whitespace-nowrap shrink-0"
              >
                Top Axis
              </button>
              <button
                type="button"
                onClick={() => onTriggerCameraPreset('focus-selected')}
                className="px-2.5 py-1 text-[11px] font-medium text-cyan-300 bg-slate-800 hover:bg-slate-700 rounded transition-colors whitespace-nowrap shrink-0"
              >
                Focus Residue
              </button>
              <button
                type="button"
                onClick={onToggleLightBackground}
                className="px-2.5 py-1 text-[11px] font-medium text-slate-300 bg-slate-800 hover:bg-slate-700 rounded transition-colors whitespace-nowrap shrink-0"
              >
                {lightBackground ? 'Dark Stage' : 'Light Stage'}
              </button>
            </div>
          </div>

          {/* Error Banner if live load failed */}
          {loadError && (
            <div
              role="alert"
              className="px-4 py-1.5 bg-rose-950/90 border-b border-rose-800 text-xs text-rose-200 flex items-center justify-between"
            >
              <span>Notice: {loadError}</span>
            </div>
          )}

          {/* TAB 1: 3D VIEWPORT & H-BOND CONTROLS */}
          {activeTab === 'viewport' && (
            <div className="p-3 sm:p-4 grid grid-cols-1 md:grid-cols-12 gap-3 sm:gap-4 items-center">
              {/* Representation & Color Scheme */}
              <div className="md:col-span-5 space-y-2.5">
                <div>
                  <span className="block text-[11px] font-medium text-slate-400 mb-1">
                    3D Structural Representation
                  </span>
                  <div className="grid grid-cols-2 sm:grid-cols-4 gap-1 bg-slate-900 p-1 rounded-lg border border-slate-800">
                    {(
                      [
                        { id: 'cartoon-sidechain', label: 'Cartoon' },
                        { id: 'ball-and-stick', label: 'Ball & Stick' },
                        { id: 'backbone-trace', label: 'Backbone' },
                        { id: 'spacefill', label: 'Spacefill' },
                      ] as const
                    ).map((item) => (
                      <button
                        key={item.id}
                        type="button"
                        onClick={() => onChangeRepresentation(item.id)}
                        className={`px-2 py-1 text-[11px] font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                          representation === item.id
                            ? 'bg-cyan-500 text-slate-950 font-semibold'
                            : 'text-slate-300 hover:text-white'
                        }`}
                      >
                        {item.label}
                      </button>
                    ))}
                  </div>
                </div>

                <div>
                  <span className="block text-[11px] font-medium text-slate-400 mb-1">
                    Color Encoding Scheme
                  </span>
                  <div className="flex flex-wrap gap-1 bg-slate-900 p-1 rounded-lg border border-slate-800">
                    {(
                      [
                        { id: 'plddt-bfactor', label: 'pLDDT / B-Factor' },
                        { id: 'secondary-structure', label: 'Secondary Struct' },
                        { id: 'element', label: 'CPK Element' },
                        { id: 'hydrophobicity', label: 'Hydrophobicity' },
                        { id: 'chain', label: 'Chain ID' },
                      ] as const
                    ).map((item) => (
                      <button
                        key={item.id}
                        type="button"
                        onClick={() => onChangeColorScheme(item.id)}
                        className={`px-2.5 py-1 text-[11px] font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                          colorScheme === item.id
                            ? 'bg-cyan-500 text-slate-950 font-semibold'
                            : 'text-slate-300 hover:text-white'
                        }`}
                      >
                        {item.label}
                      </button>
                    ))}
                  </div>
                </div>
              </div>

              {/* Hydrogen Bond Network Real-Time Controls */}
              <div className="md:col-span-4 space-y-2.5 border-t md:border-t-0 md:border-l border-slate-800 md:pl-4 pt-2 md:pt-0">
                <div className="flex items-center justify-between gap-2">
                  <label
                    htmlFor="hbond-cutoff-slider"
                    className="text-[11px] font-medium text-slate-300"
                  >
                    H-Bond Cutoff:{' '}
                    <span className="font-mono font-bold text-cyan-400 tabular-nums">
                      {hbondCutoff.toFixed(2)} Å
                    </span>{' '}
                    <span className="font-mono text-slate-400 tabular-nums">
                      ({filteredHydrogenBonds.length} bonds)
                    </span>
                  </label>
                  <button
                    type="button"
                    onClick={onToggleHydrogenBonds}
                    aria-pressed={showHydrogenBonds}
                    className={`px-2.5 py-0.5 text-[11px] font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                      showHydrogenBonds
                        ? 'bg-emerald-500/20 text-emerald-300 border border-emerald-500/50'
                        : 'bg-slate-800 text-slate-400 border border-slate-700'
                    }`}
                  >
                    {showHydrogenBonds ? 'H-Bonds Visible' : 'H-Bonds Hidden'}
                  </button>
                </div>

                <input
                  id="hbond-cutoff-slider"
                  type="range"
                  min={2.4}
                  max={4.0}
                  step={0.05}
                  value={hbondCutoff}
                  onChange={(e) => onChangeHBondCutoff(parseFloat(e.target.value))}
                  className="w-full accent-cyan-400 cursor-pointer h-1.5 bg-slate-800 rounded-lg"
                />

                <div className="flex flex-wrap gap-1 bg-slate-900 p-1 rounded-lg border border-slate-800">
                  {(
                    [
                      { id: 'all', label: 'All H-Bonds' },
                      { id: 'backbone', label: 'Backbone α/β' },
                      { id: 'sidechain-ligand', label: 'Sidechain & Ligand' },
                      { id: 'selected-only', label: 'Selected Res' },
                    ] as const
                  ).map((f) => (
                    <button
                      key={f.id}
                      type="button"
                      onClick={() => onChangeHBondFilter(f.id)}
                      className={`px-2 py-1 text-[11px] font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                        hbondFilter === f.id
                          ? 'bg-cyan-500 text-slate-950 font-semibold'
                          : 'text-slate-300 hover:text-white'
                      }`}
                    >
                      {f.label}
                    </button>
                  ))}
                </div>
              </div>

              {/* Exploded Domain Slider & Canonical AlphaFold pLDDT Scale Legend */}
              <div className="md:col-span-3 space-y-2.5 border-t md:border-t-0 md:border-l border-slate-800 md:pl-4 pt-2 md:pt-0">
                <div>
                  <div className="flex items-center justify-between text-[11px] text-slate-300 mb-1">
                    <label htmlFor="exploded-domain-slider">Domain Expansion</label>
                    <span className="font-mono tabular-nums text-cyan-400">
                      {Math.round(explodedScale * 100)}%
                    </span>
                  </div>
                  <input
                    id="exploded-domain-slider"
                    type="range"
                    min={0}
                    max={1}
                    step={0.05}
                    value={explodedScale}
                    onChange={(e) => onChangeExplodedScale(parseFloat(e.target.value))}
                    className="w-full accent-cyan-400 cursor-pointer h-1.5 bg-slate-800 rounded-lg"
                  />
                </div>

                {/* Canonical 4-Band AlphaFold pLDDT Confidence Scale */}
                <div className="pt-1">
                  <div className="text-[10px] text-slate-400 mb-1">
                    AlphaFold pLDDT Confidence Scale
                  </div>
                  <div className="grid grid-cols-2 gap-x-2 gap-y-1 text-[10px] font-mono text-slate-300">
                    <span className="flex items-center gap-1.5 whitespace-nowrap">
                      <span className="w-2.5 h-2.5 rounded-xs bg-[#0053D6] inline-block shrink-0" />
                      &gt;90 Very High
                    </span>
                    <span className="flex items-center gap-1.5 whitespace-nowrap">
                      <span className="w-2.5 h-2.5 rounded-xs bg-[#65CBF3] inline-block shrink-0" />
                      70–90 Confident
                    </span>
                    <span className="flex items-center gap-1.5 whitespace-nowrap">
                      <span className="w-2.5 h-2.5 rounded-xs bg-[#FFDB13] inline-block shrink-0" />
                      50–70 Low
                    </span>
                    <span className="flex items-center gap-1.5 whitespace-nowrap">
                      <span className="w-2.5 h-2.5 rounded-xs bg-[#FF7D45] inline-block shrink-0" />
                      &lt;50 Disordered
                    </span>
                  </div>
                </div>
              </div>
            </div>
          )}

          {/* TAB 2: INTERACTIVE RESIDUE SEQUENCE & ATOMIC COORDINATE INSPECTOR */}
          {activeTab === 'sequence' && (
            <div className="p-3.5 space-y-3">
              <div className="flex flex-wrap items-center justify-between gap-2">
                <div className="text-xs text-slate-300">
                  Click any residue in the linear sequence to highlight its sidechain rotamer,
                  inspect atomic coordinates, and focus the 3D camera.
                </div>
                <div className="flex items-center gap-2">
                  <label htmlFor="sequence-filter-input" className="sr-only">
                    Filter residues by name or number
                  </label>
                  <input
                    id="sequence-filter-input"
                    type="text"
                    value={seqFilter}
                    onChange={(e) => setSeqFilter(e.target.value)}
                    placeholder="Filter e.g. CYS, 108, HELIX..."
                    className="px-2.5 py-1 text-xs font-mono bg-slate-900 border border-slate-700 rounded text-slate-100 placeholder:text-slate-500 focus:outline-none focus:border-cyan-400"
                  />
                </div>
              </div>

              {/* Horizontal Scrollable Residue Strip */}
              <div
                role="listbox"
                aria-label="Polypeptide residue sequence"
                className="flex items-center gap-1 overflow-x-auto pb-2 pt-0.5"
              >
                {filteredResidues.map((res) => {
                  const isSelected = res.key === selectedResidue?.key;
                  const plddtColor = getPlddtColorHex(res.plddtOrBfactor);
                  const ssGlyph =
                    res.secondaryStructure === 'HELIX'
                      ? 'α'
                      : res.secondaryStructure === 'SHEET'
                      ? 'β'
                      : '—';

                  return (
                    <button
                      key={res.key}
                      role="option"
                      aria-selected={isSelected}
                      type="button"
                      onClick={() => onSelectResidue(res.key)}
                      className={`flex flex-col items-center px-2 py-1 rounded border transition-colors shrink-0 font-mono tabular-nums ${
                        isSelected
                          ? 'bg-cyan-500/20 border-cyan-400 text-white'
                          : 'bg-slate-900/90 border-slate-800 text-slate-300 hover:border-slate-600'
                      }`}
                    >
                      <span className="text-[10px] text-slate-400">{res.resSeq}</span>
                      <span className="text-xs font-bold">{res.oneLetter}</span>
                      <span className="text-[10px] text-slate-400">{res.resName}</span>
                      <span className="text-[10px] text-cyan-300">{ssGlyph}</span>
                      <span
                        className="w-full h-1 rounded-full mt-1"
                        style={{ backgroundColor: plddtColor }}
                      />
                    </button>
                  );
                })}
              </div>

              {/* Selected Residue Heavy-Atom Coordinate & Bond Summary */}
              {selectedResidue && (
                <div className="grid grid-cols-1 md:grid-cols-2 gap-3 pt-1 border-t border-slate-800 text-xs">
                  <div className="bg-slate-900/80 border border-slate-800 rounded p-2.5">
                    <div className="flex items-center justify-between mb-1.5">
                      <span className="font-semibold text-slate-200 font-mono">
                        {selectedResidue.resName} {selectedResidue.resSeq} (Chain{' '}
                        {selectedResidue.chainId}) — Heavy Atoms ({selectedResidue.atoms.length})
                      </span>
                      <span className="font-mono text-[11px] text-cyan-400 tabular-nums">
                        {model.source === 'ALPHAFOLD_DB' ? 'pLDDT' : 'B-Factor'}:{' '}
                        {selectedResidue.plddtOrBfactor.toFixed(1)}
                      </span>
                    </div>
                    <div className="max-h-24 overflow-y-auto font-mono text-[11px] tabular-nums space-y-1">
                      {selectedResidue.atoms.map((a) => (
                        <div
                          key={a.serial}
                          className="flex items-center justify-between text-slate-300 border-b border-slate-800/60 pb-0.5"
                        >
                          <span>
                            #{a.serial} <strong>{a.name}</strong> ({a.element})
                          </span>
                          <span className="text-slate-400">
                            X: {a.x.toFixed(2)} · Y: {a.y.toFixed(2)} · Z: {a.z.toFixed(2)} Å
                          </span>
                        </div>
                      ))}
                    </div>
                  </div>

                  <div className="bg-slate-900/80 border border-slate-800 rounded p-2.5">
                    <div className="font-semibold text-slate-200 mb-1.5">
                      Participating Hydrogen Bonds & Functional Context
                    </div>
                    <div className="max-h-24 overflow-y-auto space-y-1 font-mono text-[11px] tabular-nums">
                      {filteredHydrogenBonds
                        .filter(
                          (hb) =>
                            hb.donorResKey === selectedResidue.key ||
                            hb.acceptorResKey === selectedResidue.key
                        )
                        .map((hb) => (
                          <button
                            key={hb.id}
                            type="button"
                            onClick={() => onSelectHBond(hb)}
                            className="w-full text-left flex items-center justify-between px-2 py-1 rounded bg-slate-800/80 hover:bg-slate-700 text-slate-200 transition-colors"
                          >
                            <span>
                              {hb.donorAtom.resName}
                              {hb.donorAtom.resSeq}.{hb.donorAtom.name} →{' '}
                              {hb.acceptorAtom.resName}
                              {hb.acceptorAtom.resSeq}.{hb.acceptorAtom.name}
                            </span>
                            <span className="text-cyan-300">
                              {hb.distance.toFixed(2)} Å · {hb.energyKcalMol} kcal/mol
                            </span>
                          </button>
                        ))}
                      {filteredHydrogenBonds.every(
                        (hb) =>
                          hb.donorResKey !== selectedResidue.key &&
                          hb.acceptorResKey !== selectedResidue.key
                      ) && (
                        <div className="text-slate-400 py-2">
                          No hydrogen bonds connected to {selectedResidue.resName}{' '}
                          {selectedResidue.resSeq} within {hbondCutoff.toFixed(2)} Å cutoff.
                        </div>
                      )}
                    </div>
                  </div>
                </div>
              )}
            </div>
          )}

          {/* TAB 3: GEOMETRIC HYDROGEN BOND NETWORK TABLE */}
          {activeTab === 'hbonds' && (
            <div className="p-3.5 space-y-2.5">
              <div className="flex flex-wrap items-center justify-between gap-2 text-xs">
                <div className="text-slate-300">
                  Showing{' '}
                  <span className="font-mono font-bold text-cyan-400 tabular-nums">
                    {filteredHydrogenBonds.length}
                  </span>{' '}
                  donor–acceptor hydrogen bonds within{' '}
                  <span className="font-mono text-slate-100 tabular-nums">
                    {hbondCutoff.toFixed(2)} Å
                  </span>
                  . Click any bond row to inspect and center in 3D.
                </div>
                <div className="flex items-center gap-3 text-[11px] font-mono text-slate-300">
                  <span>● Cyan: α-Helix (i→i+4)</span>
                  <span>● Violet: β-Sheet Cross-Strand</span>
                  <span>● Emerald: Sidechain</span>
                  <span>● Amber: Ligand Pocket</span>
                </div>
              </div>

              <div className="max-h-44 overflow-y-auto border border-slate-800 rounded-lg">
                <table className="w-full text-left border-collapse text-xs font-mono tabular-nums">
                  <thead className="bg-slate-900 text-slate-400 text-[11px] sticky top-0">
                    <tr>
                      <th className="py-1.5 px-3">ID</th>
                      <th className="py-1.5 px-3">Donor (D–H)</th>
                      <th className="py-1.5 px-3">Acceptor (:A)</th>
                      <th className="py-1.5 px-3">Distance</th>
                      <th className="py-1.5 px-3">Angle</th>
                      <th className="py-1.5 px-3">ΔE (kcal/mol)</th>
                      <th className="py-1.5 px-3">Classification</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-slate-800/70">
                    {filteredHydrogenBonds.map((hb) => {
                      const isActive = hb.id === activeHBondId;
                      return (
                        <tr
                          key={hb.id}
                          onClick={() => onSelectHBond(hb)}
                          className={`cursor-pointer transition-colors ${
                            isActive
                              ? 'bg-cyan-500/20 text-white'
                              : 'hover:bg-slate-800/70 text-slate-200'
                          }`}
                        >
                          <td className="py-1.5 px-3 text-cyan-400 font-semibold">{hb.id}</td>
                          <td className="py-1.5 px-3">
                            {hb.donorAtom.resName} {hb.donorAtom.resSeq} · {hb.donorAtom.name} (
                            {hb.donorAtom.chainId})
                          </td>
                          <td className="py-1.5 px-3">
                            {hb.acceptorAtom.resName} {hb.acceptorAtom.resSeq} ·{' '}
                            {hb.acceptorAtom.name} ({hb.acceptorAtom.chainId})
                          </td>
                          <td className="py-1.5 px-3 font-semibold text-emerald-300">
                            {hb.distance.toFixed(2)} Å
                          </td>
                          <td className="py-1.5 px-3 text-slate-300">{hb.angleDegrees}°</td>
                          <td className="py-1.5 px-3 text-slate-300">{hb.energyKcalMol}</td>
                          <td className="py-1.5 px-3 text-slate-300">{hb.category}</td>
                        </tr>
                      );
                    })}
                    {filteredHydrogenBonds.length === 0 && (
                      <tr>
                        <td colSpan={7} className="py-6 text-center text-slate-400">
                          No hydrogen bonds match the current filter / cutoff ({hbondCutoff.toFixed(2)} Å).
                          Increase the distance cutoff slider or switch filter to "All H-Bonds".
                        </td>
                      </tr>
                    )}
                  </tbody>
                </table>
              </div>
            </div>
          )}

          {/* TAB 4: ALPHAFOLD DB pLDDT & PAE CONFIDENCE SYNTHESIS */}
          {activeTab === 'confidence' && (
            <div className="p-3.5 grid grid-cols-1 md:grid-cols-12 gap-4 text-xs">
              {/* Cohesive Structural & PAE Synthesis + Explicit Disorder Warning */}
              <div className="md:col-span-6 space-y-2">
                <div className="font-semibold text-slate-100">
                  Structural Confidence & Domain Boundary Synthesis
                </div>
                <p className="text-slate-300 leading-relaxed text-[11px]">
                  {model.alphaFoldAnalysis.cohesiveSummary}
                </p>

                {model.alphaFoldAnalysis.disorderWarning && (
                  <div
                    role="note"
                    className="p-2 rounded bg-amber-950/80 border border-amber-700/80 text-amber-200 text-[11px] leading-snug"
                  >
                    ▲ {model.alphaFoldAnalysis.disorderWarning}
                  </div>
                )}

                {model.alphaFoldAnalysis.isoformOrLengthWarning && (
                  <div className="p-2 rounded bg-slate-900 border border-slate-700 text-slate-300 text-[11px]">
                    {model.alphaFoldAnalysis.isoformOrLengthWarning}
                  </div>
                )}

                {/* 4-Band Fraction Breakdown */}
                <div className="grid grid-cols-4 gap-2 pt-1 font-mono tabular-nums">
                  <div className="p-2 rounded bg-slate-900 border border-slate-800">
                    <div className="text-[10px] text-slate-400">&gt;90 Very High</div>
                    <div className="text-sm font-bold text-[#60a5fa]">
                      {model.alphaFoldAnalysis.fractionVeryHigh}%
                    </div>
                  </div>
                  <div className="p-2 rounded bg-slate-900 border border-slate-800">
                    <div className="text-[10px] text-slate-400">70–90 Confident</div>
                    <div className="text-sm font-bold text-[#65CBF3]">
                      {model.alphaFoldAnalysis.fractionConfident}%
                    </div>
                  </div>
                  <div className="p-2 rounded bg-slate-900 border border-slate-800">
                    <div className="text-[10px] text-slate-400">50–70 Low</div>
                    <div className="text-sm font-bold text-[#FFDB13]">
                      {model.alphaFoldAnalysis.fractionLow}%
                    </div>
                  </div>
                  <div className="p-2 rounded bg-slate-900 border border-slate-800">
                    <div className="text-[10px] text-slate-400">&lt;50 Very Low</div>
                    <div className="text-sm font-bold text-[#FF7D45]">
                      {model.alphaFoldAnalysis.fractionVeryLow}%
                    </div>
                  </div>
                </div>
              </div>

              {/* PAE 2D Matrix Canvas */}
              <div className="md:col-span-3 flex flex-col items-center justify-center bg-slate-900/80 border border-slate-800 rounded p-2">
                <div className="text-[11px] font-medium text-slate-300 mb-1">
                  PAE Domain Error Matrix (Å)
                </div>
                <canvas
                  ref={paeCanvasRef}
                  width={120}
                  height={120}
                  className="border border-slate-700 rounded"
                  aria-label="Predicted Aligned Error 2D heatmap matrix"
                />
                <div className="mt-1 text-[10px] font-mono text-slate-400 flex items-center justify-between w-full px-2">
                  <span className="text-teal-400">0 Å (Rigid)</span>
                  <span className="text-amber-400">30 Å (Flexible)</span>
                </div>
              </div>

              {/* Functional Sites & Clinical Variants */}
              <div className="md:col-span-3 space-y-1.5">
                <div className="text-[11px] font-semibold text-slate-200">
                  Annotated Sites & Variants (Click to 3D Focus)
                </div>
                <div className="max-h-36 overflow-y-auto space-y-1.5">
                  {model.annotations.map((ann, idx) => (
                    <button
                      key={idx}
                      type="button"
                      onClick={() => {
                        const target = model.residues.find(
                          (r) => r.resSeq >= ann.startRes && r.resSeq <= ann.endRes
                        );
                        if (target) onSelectResidue(target.key);
                      }}
                      className="w-full text-left p-2 rounded bg-slate-900 hover:bg-slate-800 border border-slate-800 transition-colors"
                    >
                      <div className="flex items-center justify-between font-mono text-[11px] text-cyan-300">
                        <span>{ann.type}</span>
                        <span>
                          Res {ann.startRes}
                          {ann.endRes !== ann.startRes ? `–${ann.endRes}` : ''}
                        </span>
                      </div>
                      <div className="text-[11px] text-slate-300 mt-0.5 line-clamp-2">
                        {ann.description}
                      </div>
                    </button>
                  ))}
                </div>
              </div>
            </div>
          )}

          {/* TAB 5: LIVE RCSB PDB & ALPHAFOLD DB SEARCH */}
          {activeTab === 'search' && (
            <div className="p-3.5 space-y-3 text-xs">
              <form
                onSubmit={onSubmitSearch}
                className="flex flex-wrap items-center gap-2"
              >
                <label htmlFor="pdb-af-search-input" className="sr-only">
                  Search RCSB PDB ID, UniProt ID, or protein title
                </label>
                <input
                  id="pdb-af-search-input"
                  type="text"
                  value={searchQuery}
                  onChange={(e) => onSearchQueryChange(e.target.value)}
                  placeholder="Enter PDB ID (1CRN, 4HHB), UniProt ID (P04637, P01308), or keyword..."
                  className="flex-1 min-w-[240px] px-3 py-1.5 text-xs font-mono bg-slate-900 border border-slate-700 rounded-md text-slate-100 placeholder:text-slate-500 focus:outline-none focus:border-cyan-400"
                />

                <div className="flex items-center gap-1.5 text-[11px] text-slate-300 font-mono">
                  <label htmlFor="max-res-select">Max Res:</label>
                  <select
                    id="max-res-select"
                    value={maxResolutionFilter}
                    onChange={(e) => onMaxResolutionChange(parseFloat(e.target.value))}
                    className="bg-slate-900 border border-slate-700 rounded px-2 py-1 text-slate-200"
                  >
                    <option value={1.5}>≤ 1.5 Å</option>
                    <option value={2.0}>≤ 2.0 Å</option>
                    <option value={2.5}>≤ 2.5 Å</option>
                    <option value={3.5}>≤ 3.5 Å</option>
                  </select>
                </div>

                <button
                  type="submit"
                  disabled={isSearching}
                  className="px-3.5 py-1.5 text-xs font-semibold text-slate-950 bg-cyan-400 hover:bg-cyan-300 rounded-md transition-colors whitespace-nowrap shrink-0"
                >
                  {isSearching ? 'Querying APIs...' : 'Search Databases'}
                </button>

                <button
                  type="button"
                  disabled={isLoadingStructure || !searchQuery.trim()}
                  onClick={() => onLoadAccession(searchQuery.trim(), true)}
                  className="px-3 py-1.5 text-xs font-medium text-cyan-300 bg-slate-800 hover:bg-slate-700 rounded-md transition-colors whitespace-nowrap shrink-0"
                >
                  Direct Fetch Accession
                </button>
              </form>

              {searchError && (
                <div role="alert" className="text-rose-400 text-xs font-medium">
                  ✖ {searchError}
                </div>
              )}

              <div className="text-[11px] text-slate-400 font-mono">
                {searchExplanation || queryExplanation}
              </div>

              <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 gap-2 max-h-36 overflow-y-auto">
                {searchResults.map((item) => (
                  <button
                    key={`${item.source}-${item.id}`}
                    type="button"
                    onClick={() => onLoadAccession(item.id)}
                    className="text-left p-2.5 rounded-lg bg-slate-900/90 hover:bg-slate-800 border border-slate-800 transition-colors flex flex-col justify-between"
                  >
                    <div className="flex items-center justify-between gap-2 font-mono text-xs">
                      <span className="font-bold text-cyan-400">{item.id}</span>
                      <span className="text-[10px] text-slate-400">
                        {item.source === 'ALPHAFOLD_DB' ? 'AlphaFold DB' : 'RCSB PDB'}
                        {item.resolution ? ` · ${item.resolution.toFixed(2)} Å` : ''}
                      </span>
                    </div>
                    <div className="text-xs font-medium text-slate-200 truncate mt-0.5">
                      {item.title}
                    </div>
                    <div className="text-[10px] text-slate-400 truncate mt-0.5">
                      {item.method} · {item.organism}
                    </div>
                  </button>
                ))}
              </div>
            </div>
          )}
        </div>
      )}
    </section>
  );
};
