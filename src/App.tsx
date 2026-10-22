/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useCallback, useEffect, useMemo, useState } from 'react';
import { BottomControlDock } from './components/BottomControlDock';
import { CameraPresetType, ProteinViewport3D } from './components/ProteinViewport3D';
import { DockTabId, TopBar } from './components/TopBar';
import {
   loadProteinStructure,
  searchStructuresCatalog,
} from './services/proteinApi';
import {
  AtomRecord,
  ColorSchemeMode,
  HBondFilterMode,
  HydrogenBond,
  PdbSearchResultItem,
  ProteinStructureModel,
  RepresentationMode,
} from './types/protein';
import { SEED_PROTEIN_MODELS } from './utils/seedProteins';

const INITIAL_MODEL = SEED_PROTEIN_MODELS['P04637'];

const INITIAL_SEARCH_RESULTS: PdbSearchResultItem[] = Object.values(SEED_PROTEIN_MODELS)
  .filter((m, idx, arr) => arr.findIndex((x) => x.displayCode === m.displayCode) === idx)
  .map((m) => ({
    id: m.displayCode,
    source: m.source,
    title: m.title,
    method: m.experimentalMethod,
    resolution: m.resolutionAngstrom,
    organism: m.organism,
    gene: m.geneSymbol,
  }));

export default function App() {
  const [model, setModel] = useState<ProteinStructureModel>(INITIAL_MODEL);
  const [activeTab, setActiveTab] = useState<DockTabId>('viewport');
  const [dockCollapsed, setDockCollapsed] = useState<boolean>(false);

  const [representation, setRepresentation] =
    useState<RepresentationMode>('cartoon-sidechain');
  const [colorScheme, setColorScheme] = useState<ColorSchemeMode>('plddt-bfactor');

  const [showHydrogenBonds, setShowHydrogenBonds] = useState<boolean>(true);
  const [hbondCutoff, setHbondCutoff] = useState<number>(3.35);
  const [hbondFilter, setHbondFilter] = useState<HBondFilterMode>('all');

  const [selectedResidueKey, setSelectedResidueKey] = useState<string | null>(
    INITIAL_MODEL.residues[12]?.key || INITIAL_MODEL.residues[0]?.key || null
  );
  const [hoveredResidueKey, setHoveredResidueKey] = useState<string | null>(null);
  const [hoveredAtom, setHoveredAtom] = useState<AtomRecord | null>(null);
  const [activeHBondId, setActiveHBondId] = useState<string | null>(null);

  const [lightBackground, setLightBackground] = useState<boolean>(false);
  const [explodedScale, setExplodedScale] = useState<number>(0);
  const [cameraPresetTrigger, setCameraPresetTrigger] = useState<{
    preset: CameraPresetType;
    seq: number;
  }>({ preset: 'overview', seq: 1 });

  const [queryExplanation, setQueryExplanation] = useState<string>(
    'Resolved AlphaFold DB prediction metadata via GET https://alphafold.ebi.ac.uk/api/prediction/P04637 (dynamic _v6 coordinate URL) and computed per-residue pLDDT confidence from the B-factor column.'
  );
  const [searchQuery, setSearchQuery] = useState<string>('TP53');
  const [maxResolutionFilter, setMaxResolutionFilter] = useState<number>(2.5);
  const [searchResults, setSearchResults] =
    useState<PdbSearchResultItem[]>(INITIAL_SEARCH_RESULTS);
  const [searchExplanation, setSearchExplanation] = useState<string>('');
  const [searchError, setSearchError] = useState<string | null>(null);
  const [isSearching, setIsSearching] = useState<boolean>(false);

  const [isLoadingStructure, setIsLoadingStructure] = useState<boolean>(false);
  const [loadError, setLoadError] = useState<string | null>(null);
  const [exportNotice, setExportNotice] = useState<string | null>(null);

  // Keyboard shortcut 'H' to toggle collapsible Zen Mode Bottom Dock
  useEffect(() => {
    const onKeyDown = (e: KeyboardEvent) => {
      const tag = (e.target as HTMLElement)?.tagName;
      if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return;
      if (e.key === 'h' || e.key === 'H') {
        setDockCollapsed((prev) => !prev);
      }
    };
    window.addEventListener('keydown', onKeyDown);
    return () => window.removeEventListener('keydown', onKeyDown);
  }, []);

  // Filter Hydrogen Bonds in real-time by cutoff distance and category
  const filteredHydrogenBonds = useMemo(() => {
    return model.allHydrogenBonds.filter((hb) => {
      if (hb.distance > hbondCutoff) return false;
      if (hbondFilter === 'backbone') {
        return hb.category === 'backbone-helix' || hb.category === 'backbone-sheet';
      }
      if (hbondFilter === 'sidechain-ligand') {
        return hb.category === 'sidechain' || hb.category === 'ligand';
      }
      if (hbondFilter === 'selected-only' && selectedResidueKey) {
        return (
          hb.donorResKey === selectedResidueKey ||
          hb.acceptorResKey === selectedResidueKey
        );
      }
      return true;
    });
  }, [model.allHydrogenBonds, hbondCutoff, hbondFilter, selectedResidueKey]);

  const selectedResidue = useMemo(
    () => model.residues.find((r) => r.key === selectedResidueKey) || null,
    [model.residues, selectedResidueKey]
  );

  const hoveredResidue = useMemo(
    () => model.residues.find((r) => r.key === hoveredResidueKey) || null,
    [model.residues, hoveredResidueKey]
  );

  const selectedResidueHBonds = useMemo(() => {
    const targetKey = hoveredResidueKey || selectedResidueKey;
    if (!targetKey) return [];
    return filteredHydrogenBonds.filter(
      (hb) => hb.donorResKey === targetKey || hb.acceptorResKey === targetKey
    );
  }, [filteredHydrogenBonds, hoveredResidueKey, selectedResidueKey]);

  const triggerCameraPreset = useCallback((preset: CameraPresetType) => {
    setCameraPresetTrigger((prev) => ({ preset, seq: prev.seq + 1 }));
  }, []);

  const handleSelectResidue = useCallback(
    (residueKey: string) => {
      setSelectedResidueKey(residueKey);
      setActiveHBondId(null);
      triggerCameraPreset('focus-selected');
    },
    [triggerCameraPreset]
  );

  const handleHoverResidue = useCallback(
    (residueKey: string | null, atomInfo?: AtomRecord | null) => {
      setHoveredResidueKey(residueKey);
      setHoveredAtom(atomInfo || null);
    },
    []
  );

  const handleSelectHBond = useCallback(
    (hbond: HydrogenBond) => {
      setActiveHBondId(hbond.id);
      setSelectedResidueKey(hbond.donorResKey);
      setShowHydrogenBonds(true);
      triggerCameraPreset('focus-hbonds');
    },
    [triggerCameraPreset]
  );

  const handleLoadAccession = useCallback(
    async (accessionCode: string, forceLive = false) => {
      if (!accessionCode.trim()) {
        setLoadError('Please specify a valid PDB ID (e.g. 1CRN) or UniProt ID (e.g. P04637).');
        return;
      }
      setIsLoadingStructure(true);
      setLoadError(null);
      try {
        const { model: loadedModel, queryExplanation: explanation } =
          await loadProteinStructure(accessionCode, forceLive);
        setModel(loadedModel);
        setQueryExplanation(explanation);
        const defaultRes =
          loadedModel.residues[Math.min(10, loadedModel.residues.length - 1)] ||
          loadedModel.residues[0];
        setSelectedResidueKey(defaultRes ? defaultRes.key : null);
        setActiveHBondId(null);
        setExplodedScale(0);
        triggerCameraPreset('overview');
      } catch (err: any) {
        setLoadError(
          err?.message || `Failed to load structure "${accessionCode}".`
        );
      } finally {
        setIsLoadingStructure(false);
      }
    },
    [triggerCameraPreset]
  );

  const handleSubmitSearch = useCallback(
    async (e: React.FormEvent) => {
      e.preventDefault();
      setSearchError(null);
      if (!searchQuery.trim()) {
        setSearchError(
          'Please enter a protein name, gene symbol, PDB ID (e.g. 1CRN), or UniProt ID (e.g. P04637).'
        );
        return;
      }
      setIsSearching(true);
      try {
        const { results, explanation } = await searchStructuresCatalog(
          searchQuery,
          maxResolutionFilter
        );
        setSearchResults(results);
        setSearchExplanation(explanation);
      } catch (err: any) {
        setSearchError(err?.message || 'Database search failed.');
      } finally {
        setIsSearching(false);
      }
    },
    [searchQuery, maxResolutionFilter]
  );

  const handleExportReport = useCallback(() => {
    const payload = {
      exportedAt: new Date().toISOString(),
      structureId: model.id,
      title: model.title,
      source: model.source,
      provenance: model.provenanceBadge,
      experimentalMethod: model.experimentalMethod,
      resolutionAngstrom: model.resolutionAngstrom,
      residueCount: model.residues.length,
      atomCount: model.atoms.length,
      hbondDistanceCutoffAngstrom: hbondCutoff,
      activeHydrogenBondCount: filteredHydrogenBonds.length,
      alphaFoldAnalysis: model.alphaFoldAnalysis,
      selectedResidue: selectedResidue
        ? {
            key: selectedResidue.key,
            resName: selectedResidue.resName,
            resSeq: selectedResidue.resSeq,
            secondaryStructure: selectedResidue.secondaryStructure,
            plddtOrBfactor: selectedResidue.plddtOrBfactor,
          }
        : null,
      hydrogenBonds: filteredHydrogenBonds.map((hb) => ({
        id: hb.id,
        donor: `${hb.donorAtom.resName}${hb.donorAtom.resSeq}.${hb.donorAtom.name}`,
        acceptor: `${hb.acceptorAtom.resName}${hb.acceptorAtom.resSeq}.${hb.acceptorAtom.name}`,
        distanceAngstrom: hb.distance,
        angleDegrees: hb.angleDegrees,
        energyKcalMol: hb.energyKcalMol,
        category: hb.category,
      })),
    };

    const blob = new Blob([JSON.stringify(payload, null, 2)], {
      type: 'application/json',
    });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${model.id.toLowerCase()}-hbond-analysis.json`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);

    setExportNotice(
      `Exported ${filteredHydrogenBonds.length} hydrogen bonds & coordinates for ${model.id}.`
    );
    window.setTimeout(() => setExportNotice(null), 4000);
  }, [model, hbondCutoff, filteredHydrogenBonds, selectedResidue]);

  const handleSelectTab = useCallback((tab: DockTabId) => {
    setActiveTab(tab);
    setDockCollapsed(false);
  }, []);

  return (
    <main className="relative w-screen h-screen overflow-hidden select-none bg-[#0b0f19]">
      {/* Full-Bleed 3D Three.js Macromolecular Viewport */}
      <ProteinViewport3D
        model={model}
        representation={representation}
        colorScheme={colorScheme}
        showHydrogenBonds={showHydrogenBonds}
        filteredHydrogenBonds={filteredHydrogenBonds}
        selectedResidueKey={selectedResidueKey}
        hoveredResidueKey={hoveredResidueKey}
        activeHBondId={activeHBondId}
        lightBackground={lightBackground}
        explodedScale={explodedScale}
        cameraPresetTrigger={cameraPresetTrigger}
        onSelectResidue={handleSelectResidue}
        onHoverResidue={handleHoverResidue}
      />

      {/* Top 3-Zone Header + Active Structure Spec & Selected Residue Readout */}
      <TopBar
        model={model}
        activeTab={activeTab}
        dockCollapsed={dockCollapsed}
        selectedResidue={selectedResidue}
        hoveredResidue={hoveredResidue}
        hoveredAtom={hoveredAtom}
        selectedResidueHBonds={selectedResidueHBonds}
        filteredHBondCount={filteredHydrogenBonds.length}
        isLoading={isLoadingStructure}
        onSelectTab={handleSelectTab}
        onSelectPreset={(code) => handleLoadAccession(code, false)}
        onExportPdbReport={handleExportReport}
        onFocusSelectedResidue={() => triggerCameraPreset('focus-selected')}
      />

      {/* Export Confirmation Notice */}
      {exportNotice && (
        <div
          role="status"
          className="fixed top-14 left-1/2 -translate-x-1/2 z-30 bg-emerald-950/95 border border-emerald-600 text-emerald-200 px-4 py-1.5 rounded-md text-xs font-mono shadow-lg"
        >
          ● {exportNotice}
        </div>
      )}

      {/* Unified Collapsible Bottom Control Dock (Zero Viewport Center Occlusion) */}
      <BottomControlDock
        model={model}
        activeTab={activeTab}
        dockCollapsed={dockCollapsed}
        representation={representation}
        colorScheme={colorScheme}
        showHydrogenBonds={showHydrogenBonds}
        hbondCutoff={hbondCutoff}
        hbondFilter={hbondFilter}
        filteredHydrogenBonds={filteredHydrogenBonds}
        selectedResidueKey={selectedResidueKey}
        activeHBondId={activeHBondId}
        lightBackground={lightBackground}
        explodedScale={explodedScale}
        queryExplanation={queryExplanation}
        searchQuery={searchQuery}
        maxResolutionFilter={maxResolutionFilter}
        searchResults={searchResults}
        searchExplanation={searchExplanation}
        searchError={searchError}
        isSearching={isSearching}
        isLoadingStructure={isLoadingStructure}
        loadError={loadError}
        onSelectTab={handleSelectTab}
        onToggleCollapse={() => setDockCollapsed((prev) => !prev)}
        onChangeRepresentation={setRepresentation}
        onChangeColorScheme={setColorScheme}
        onToggleHydrogenBonds={() => setShowHydrogenBonds((prev) => !prev)}
        onChangeHBondCutoff={setHbondCutoff}
        onChangeHBondFilter={setHbondFilter}
        onSelectResidue={handleSelectResidue}
        onSelectHBond={handleSelectHBond}
        onToggleLightBackground={() => setLightBackground((prev) => !prev)}
        onChangeExplodedScale={setExplodedScale}
        onTriggerCameraPreset={triggerCameraPreset}
        onSearchQueryChange={setSearchQuery}
        onMaxResolutionChange={setMaxResolutionFilter}
        onSubmitSearch={handleSubmitSearch}
        onLoadAccession={handleLoadAccession}
      />
    </main>
  );
}
