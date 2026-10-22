import React from 'react';
import { AtomRecord, HydrogenBond, ProteinStructureModel, ResidueRecord } from '../types/protein';
import { getPlddtColorHex, getPlddtTierLabel } from '../utils/proteinParser';
import { QUICK_TRY_PRESETS } from '../utils/seedProteins';

export type DockTabId = 'viewport' | 'sequence' | 'hbonds' | 'confidence' | 'search';

interface TopBarProps {
  model: ProteinStructureModel;
  activeTab: DockTabId;
  dockCollapsed: boolean;
  selectedResidue: ResidueRecord | null;
  hoveredResidue: ResidueRecord | null;
  hoveredAtom: AtomRecord | null;
  selectedResidueHBonds: HydrogenBond[];
  filteredHBondCount: number;
  isLoading: boolean;
  onSelectTab: (tab: DockTabId) => void;
  onSelectPreset: (code: string) => void;
  onExportPdbReport: () => void;
  onFocusSelectedResidue: () => void;
}

export const TopBar: React.FC<TopBarProps> = ({
  model,
  activeTab,
  dockCollapsed,
  selectedResidue,
  hoveredResidue,
  hoveredAtom,
  selectedResidueHBonds,
  filteredHBondCount,
  isLoading,
  onSelectTab,
  onSelectPreset,
  onExportPdbReport,
  onFocusSelectedResidue,
}) => {
  const inspectedRes = hoveredResidue || selectedResidue;

  return (
    <div className="fixed top-0 left-0 right-0 z-20 pointer-events-none">
      {/* Strict 3-Zone Top Bar Contract */}
      <header className="pointer-events-auto flex items-center justify-between gap-8 px-6 py-3 bg-[#0b0f19]/95 border-b border-slate-800">
        {/* Zone 1: Single text element wordmark */}
        <a
          href="#top"
          onClick={(e) => {
            e.preventDefault();
            onSelectTab('viewport');
          }}
          className="text-base font-bold tracking-tight text-slate-100 whitespace-nowrap shrink-0"
        >
          Structura3D
        </a>

        {/* Zone 2: 5 concise single-line navigation links */}
        <nav
          aria-label="Workspace sections"
          className="hidden md:flex items-center gap-6 text-xs font-medium text-slate-300"
        >
          {(
            [
              { id: 'viewport', label: '3D Controls' },
              { id: 'sequence', label: 'Residue Sequence' },
              { id: 'hbonds', label: 'H-Bond Network' },
              { id: 'confidence', label: 'pLDDT & PAE' },
              { id: 'search', label: 'PDB & AlphaFold Search' },
            ] as const
          ).map((item) => (
            <button
              key={item.id}
              type="button"
              onClick={() => onSelectTab(item.id)}
              className={`py-1 transition-colors whitespace-nowrap shrink-0 border-b-2 ${
                activeTab === item.id && !dockCollapsed
                  ? 'text-cyan-400 border-cyan-400 font-semibold'
                  : 'text-slate-300 border-transparent hover:text-white'
              }`}
            >
              {item.label}
            </button>
          ))}
        </nav>

        {/* Zone 3: 1 Primary Action */}
        <div className="flex items-center gap-3 shrink-0">
          <button
            type="button"
            onClick={onExportPdbReport}
            className="px-3.5 py-1.5 text-xs font-semibold text-slate-950 bg-cyan-400 rounded-md hover:bg-cyan-300 transition-colors whitespace-nowrap shrink-0"
          >
            Export H-Bond JSON
          </button>
        </div>
      </header>

      {/* Compact Top-Left / Top-Center Spec & Quick-Try Bar (Leaves 3D Center Unobstructed) */}
      {!dockCollapsed && (
        <div className="px-3 sm:px-6 pt-2 flex flex-col sm:flex-row items-stretch sm:items-start justify-between gap-2 sm:gap-3">
          {/* Top-Left: Active Structure Provenance & Quick-Try Presets */}
          <div className="pointer-events-auto bg-[#0f172a]/95 border border-slate-800 rounded-lg px-3 py-2 sm:px-3.5 sm:py-2.5 max-w-2xl">
            <div className="flex flex-wrap items-center gap-1.5 text-[11px] text-slate-400 mb-0.5">
              <span className="font-mono font-semibold text-cyan-400">{model.id}</span>
              <span aria-hidden="true">·</span>
              <span className="text-slate-200 font-medium">{model.provenanceBadge}</span>
              <span className="hidden sm:inline" aria-hidden="true">·</span>
              <span className="hidden sm:inline">{model.experimentalMethod}</span>
              {model.resolutionAngstrom !== null && (
                <>
                  <span aria-hidden="true">·</span>
                  <span className="font-mono tabular-nums text-emerald-400">
                    {model.resolutionAngstrom.toFixed(2)} Å
                  </span>
                </>
              )}
              <span className="hidden md:inline" aria-hidden="true">·</span>
              <span className="hidden md:inline">{model.organism}</span>
              <span aria-hidden="true">·</span>
              <a
                href={model.sourceUrl}
                target="_blank"
                rel="noopener noreferrer"
                className="inline-flex items-center gap-1 text-cyan-400 hover:underline font-mono whitespace-nowrap"
              >
                <span>Source Record</span>
              </a>
            </div>

            <div className="flex items-center justify-between gap-3">
              <div>
                <h1 className="text-xs sm:text-sm font-semibold text-slate-100 leading-snug">
                  {model.title}
                </h1>
                <p className="text-[11px] text-slate-400 font-mono tabular-nums mt-0.5">
                  {model.residues.length} Res · {model.atoms.length} Atoms ·{' '}
                  <span className="text-cyan-300 font-semibold">
                    {filteredHBondCount} H-Bonds
                  </span>{' '}
                  · Mean {model.source === 'ALPHAFOLD_DB' ? 'pLDDT' : 'B'}:{' '}
                  <span className="text-slate-200 font-semibold">
                    {model.alphaFoldAnalysis.meanPlddt}
                  </span>
                </p>
              </div>
            </div>

            {/* Segmented Quick-Try Preset Switcher */}
            <div
              role="group"
              aria-label="Quick-try protein structures"
              className="mt-1.5 pt-1.5 border-t border-slate-800/80 flex items-center gap-1 overflow-x-auto"
            >
              <span className="text-[11px] text-slate-400 mr-1 whitespace-nowrap shrink-0">
                Presets:
              </span>
              {QUICK_TRY_PRESETS.map((preset) => {
                const active = model.displayCode === preset.code;
                return (
                  <button
                    key={preset.code}
                    type="button"
                    disabled={isLoading}
                    onClick={() => onSelectPreset(preset.code)}
                    className={`px-2 py-0.5 sm:px-2.5 sm:py-1 text-[11px] font-mono rounded transition-colors whitespace-nowrap shrink-0 ${
                      active
                        ? 'bg-cyan-500 text-slate-950 font-semibold'
                        : 'bg-slate-800/90 text-slate-300 hover:bg-slate-700 hover:text-white'
                    }`}
                  >
                    {preset.label}
                  </button>
                );
              })}
            </div>
          </div>

          {/* Top-Right: Compact Real-Time Residue Telemetry Readout */}
          {inspectedRes && (
            <div
              data-testid="selected-residue-readout"
              className="hidden sm:block pointer-events-auto bg-[#0f172a]/95 border border-slate-800 rounded-lg px-3.5 py-2.5 min-w-[265px] max-w-sm"
            >
              <div className="flex items-center justify-between gap-3">
                <div className="text-[11px] text-slate-400">
                  {hoveredResidue ? 'Hovering Residue' : 'Selected Residue'}
                  <span aria-hidden="true"> · </span>
                  <span className="font-mono text-slate-300">
                    Chain {inspectedRes.chainId}
                  </span>
                </div>
                <button
                  type="button"
                  onClick={onFocusSelectedResidue}
                  className="px-2 py-0.5 text-[11px] font-medium text-cyan-300 bg-slate-800 hover:bg-slate-700 rounded transition-colors whitespace-nowrap shrink-0"
                >
                  Center in 3D
                </button>
              </div>

              <div className="mt-1 flex items-baseline justify-between gap-3">
                <div className="text-base font-bold font-mono text-white tabular-nums">
                  {inspectedRes.resName} {inspectedRes.resSeq}{' '}
                  <span className="text-xs font-normal text-slate-400">
                    ({inspectedRes.oneLetter})
                  </span>
                </div>
                <div className="text-xs font-mono tabular-nums flex items-center gap-1.5">
                  <span
                    className="inline-block w-2.5 h-2.5 rounded-sm"
                    style={{
                      backgroundColor: getPlddtColorHex(inspectedRes.plddtOrBfactor),
                    }}
                  />
                  <span className="text-slate-200 font-semibold">
                    {inspectedRes.plddtOrBfactor.toFixed(1)}
                  </span>
                  <span className="text-[11px] text-slate-400">
                    {model.source === 'ALPHAFOLD_DB' ? 'pLDDT' : 'B'}
                  </span>
                </div>
              </div>

              <div className="mt-1 text-[11px] text-slate-400 flex flex-wrap items-center gap-1.5 font-mono tabular-nums">
                <span>SS: {inspectedRes.secondaryStructure}</span>
                <span aria-hidden="true">·</span>
                <span>
                  Cα ({inspectedRes.caAtom.x.toFixed(1)}, {inspectedRes.caAtom.y.toFixed(1)},{' '}
                  {inspectedRes.caAtom.z.toFixed(1)} Å)
                </span>
                {hoveredAtom && (
                  <>
                    <span aria-hidden="true">·</span>
                    <span className="text-amber-300">Atom {hoveredAtom.name}</span>
                  </>
                )}
              </div>

              <div className="mt-1.5 pt-1.5 border-t border-slate-800/80 text-[11px] text-slate-300 font-mono tabular-nums">
                {selectedResidueHBonds.length > 0 ? (
                  <span>
                    Participates in{' '}
                    <strong className="text-cyan-300">
                      {selectedResidueHBonds.length} H-bond
                      {selectedResidueHBonds.length === 1 ? '' : 's'}
                    </strong>{' '}
                    ({selectedResidueHBonds[0].donorAtom.resName}
                    {selectedResidueHBonds[0].donorAtom.resSeq}{' '}
                    {selectedResidueHBonds[0].donorAtom.name} -{' '}
                    {selectedResidueHBonds[0].acceptorAtom.resName}
                    {selectedResidueHBonds[0].acceptorAtom.resSeq}{' '}
                    {selectedResidueHBonds[0].acceptorAtom.name},{' '}
                    {selectedResidueHBonds[0].distance.toFixed(2)} Å)
                  </span>
                ) : (
                  <span className="text-slate-400">
                    0 hydrogen bonds at current distance cutoff ·{' '}
                    {getPlddtTierLabel(inspectedRes.plddtOrBfactor)}
                  </span>
                )}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
};
