export type SecondaryStructureType = 'HELIX' | 'SHEET' | 'COIL';

export type StructureSource = 'RCSB_PDB' | 'ALPHAFOLD_DB';

export type RepresentationMode =
  | 'cartoon-sidechain'
  | 'ball-and-stick'
  | 'backbone-trace'
  | 'spacefill';

export type ColorSchemeMode =
  | 'plddt-bfactor'
  | 'secondary-structure'
  | 'chain'
  | 'element'
  | 'hydrophobicity';

export type HBondCategory =
  | 'backbone-helix'
  | 'backbone-sheet'
  | 'sidechain'
  | 'ligand';

export type HBondFilterMode = 'all' | 'backbone' | 'sidechain-ligand' | 'selected-only';

export interface AtomRecord {
  serial: number;
  name: string;        // e.g. "CA", "N", "O", "CB", "OG1"
  element: string;     // e.g. "C", "N", "O", "S", "FE"
  resName: string;     // e.g. "THR", "PRO", "CYS", "HEM"
  chainId: string;     // e.g. "A", "B"
  resSeq: number;      // e.g. 1, 2, ...
  x: number;
  y: number;
  z: number;
  bFactor: number;     // In AlphaFold DB, stores per-residue pLDDT (0 - 100)
  isHet: boolean;      // true for HETATM (ligands, cofactors, ions)
}

export interface ResidueRecord {
  key: string;         // `${chainId}:${resSeq}`
  chainId: string;
  resSeq: number;
  resName: string;     // 3-letter code e.g. "ARG"
  oneLetter: string;   // 1-letter code e.g. "R"
  secondaryStructure: SecondaryStructureType;
  plddtOrBfactor: number;
  caAtom: AtomRecord;
  atoms: AtomRecord[];
  isLigand: boolean;
}

export interface CovalentBond {
  atomIndexA: number;
  atomIndexB: number;
  order: 1 | 2 | 'disulfide';
}

export interface HydrogenBond {
  id: string;
  donorAtom: AtomRecord;
  acceptorAtom: AtomRecord;
  donorResKey: string;
  acceptorResKey: string;
  distance: number;          // in Angstroms (Å)
  angleDegrees: number;      // estimated N-H...O / O-H...O angle in degrees
  energyKcalMol: number;     // estimated electrostatic stabilization energy
  category: HBondCategory;
}

export interface PaeDomainBlock {
  id: string;
  label: string;
  startRes: number;
  endRes: number;
  meanInternalPae: number;
  meanPlddt: number;
}

export interface AlphaFoldConfidenceAnalysis {
  meanPlddt: number;
  fractionVeryHigh: number;  // > 90
  fractionConfident: number; // 70 - 90
  fractionLow: number;       // 50 - 70
  fractionVeryLow: number;   // < 50
  plddtConclusion: string;
  paeConclusion: string;
  cohesiveSummary: string;
  disorderWarning: string | null;
  isoformOrLengthWarning: string | null;
  domains: PaeDomainBlock[];
  paeGridSample: number[][]; // Downsampled 24x24 PAE heatmap matrix (Å)
}

export interface FunctionalAnnotation {
  type: 'Active site' | 'Binding site' | 'Domain' | 'Natural variant' | 'Disulfide bond';
  startRes: number;
  endRes: number;
  description: string;
  clinicalSignificance?: 'Pathogenic' | 'Uncertain Significance' | 'Benign';
}

export interface ProteinStructureModel {
  id: string;                 // e.g. "1CRN", "AF-P04637-F1"
  displayCode: string;        // e.g. "1CRN", "P04637"
  title: string;
  subtitle: string;
  source: StructureSource;
  provenanceBadge: 'Experimentally determined' | 'Predicted by model';
  geneSymbol: string;
  organism: string;
  uniprotId?: string;
  pdbId?: string;
  experimentalMethod: string; // e.g. "X-RAY DIFFRACTION" or "AlphaFold Monomer v2.0 (v6)"
  resolutionAngstrom: number | null;
  releaseDate: string;
  citationTitle: string;
  sourceUrl: string;
  cifDownloadUrl: string;
  atoms: AtomRecord[];
  residues: ResidueRecord[];
  chains: string[];
  covalentBonds: CovalentBond[];
  allHydrogenBonds: HydrogenBond[]; // Precomputed up to 4.0 Å so slider filters in real time
  alphaFoldAnalysis: AlphaFoldConfidenceAnalysis;
  annotations: FunctionalAnnotation[];
}

export interface PdbSearchResultItem {
  id: string;
  source: StructureSource;
  title: string;
  method: string;
  resolution: number | null;
  organism: string;
  gene?: string;
}
