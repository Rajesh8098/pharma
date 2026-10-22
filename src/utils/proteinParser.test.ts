import test from 'node:test';
import assert from 'node:assert/strict';
import {
  getPlddtColorHex,
  parsePdbText,
  parseCifText,
  buildResiduesFromAtoms,
  computeHydrogenBonds,
  computeCovalentBonds,
  analyzeAlphaFoldConfidence,
} from './proteinParser';
import { SEED_PROTEIN_MODELS } from './seedProteins';

test('getPlddtColorHex maps canonical 4-band AlphaFold confidence colors', () => {
  assert.equal(getPlddtColorHex(95.2), '#0053D6');
  assert.equal(getPlddtColorHex(81.0), '#65CBF3');
  assert.equal(getPlddtColorHex(58.4), '#FFDB13');
  assert.equal(getPlddtColorHex(42.1), '#FF7D45');
});

test('SEED_PROTEIN_MODELS contain valid residues, covalent bonds, and hydrogen bonds', () => {
  const tp53 = SEED_PROTEIN_MODELS['P04637'];
  assert.ok(tp53.residues.length >= 50, 'TP53 has >= 50 residues');
  assert.ok(tp53.covalentBonds.length > 100, 'TP53 has covalent bonds');
  assert.ok(tp53.allHydrogenBonds.length >= 10, 'TP53 has detected hydrogen bonds');
  assert.ok(
    tp53.allHydrogenBonds.some((hb) => hb.category === 'backbone-helix'),
    'TP53 has alpha-helical backbone H-bonds'
  );
  assert.ok(
    tp53.allHydrogenBonds.some((hb) => hb.category === 'backbone-sheet'),
    'TP53 has beta-sheet cross-strand H-bonds'
  );
  assert.ok(
    tp53.alphaFoldAnalysis.disorderWarning !== null,
    'TP53 surfaces explicit disorder warning for flanking low-pLDDT segments'
  );

  const crn = SEED_PROTEIN_MODELS['1CRN'];
  assert.equal(crn.residues.length, 46, '1CRN Crambin has 46 residues');
  assert.ok(crn.allHydrogenBonds.length >= 8, '1CRN has >= 8 hydrogen bonds');

  const hhb = SEED_PROTEIN_MODELS['4HHB'];
  assert.ok(
    hhb.residues.some((r) => r.isLigand && r.resName === 'HEM'),
    '4HHB includes HEM prosthetic group'
  );
});

test('parsePdbText and computeHydrogenBonds respect distance cutoff', () => {
  const samplePdb = [
    'TITLE     TEST HELICAL SEGMENT',
    'REMARK   2 RESOLUTION.    1.50 ANGSTROMS.',
    'HELIX    1  H1 ALA A    1  ALA A    5  1                                   5',
    'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 92.00           N',
    'ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 92.00           C',
    'ATOM      3  C   ALA A   1       2.000   1.420   0.000  1.00 92.00           C',
    'ATOM      4  O   ALA A   1       3.210   1.650   0.000  1.00 92.00           O',
    'ATOM      5  N   ALA A   5       6.050   1.720   0.000  1.00 94.00           N',
    'ATOM      6  CA  ALA A   5       7.450   1.900   0.000  1.00 94.00           C',
    'ATOM      7  C   ALA A   5       8.100   3.200   0.000  1.00 94.00           C',
    'ATOM      8  O   ALA A   5       9.300   3.300   0.000  1.00 94.00           O',
  ].join('\n');

  const parsed = parsePdbText(samplePdb);
  assert.equal(parsed.atoms.length, 8);
  assert.equal(parsed.resolution, 1.5);
  const residues = buildResiduesFromAtoms(parsed.atoms, parsed.ssMap);
  assert.equal(residues.length, 2);

  const hbStrict = computeHydrogenBonds(parsed.atoms, residues, 2.7);
  const hbRelaxed = computeHydrogenBonds(parsed.atoms, residues, 3.2);
  assert.equal(hbStrict.length, 0, 'O(1)..N(5) at ~2.84 Å should be filtered out at 2.7 Å cutoff');
  assert.equal(hbRelaxed.length, 1, 'O(1)..N(5) at ~2.84 Å should be detected at 3.2 Å cutoff');
  assert.equal(hbRelaxed[0].category, 'backbone-helix');
});

test('parseCifText parses mmCIF _atom_site loop records accurately', () => {
  const sampleCif = [
    'data_TEST',
    'loop_',
    '_atom_site.group_PDB',
    '_atom_site.id',
    '_atom_site.type_symbol',
    '_atom_site.label_atom_id',
    '_atom_site.label_alt_id',
    '_atom_site.label_comp_id',
    '_atom_site.auth_asym_id',
    '_atom_site.auth_seq_id',
    '_atom_site.Cartn_x',
    '_atom_site.Cartn_y',
    '_atom_site.Cartn_z',
    '_atom_site.B_iso_or_equiv',
    '_atom_site.pdbx_PDB_model_num',
    'ATOM 1 N N . THR A 1 1.000 2.000 3.000 95.40 1',
    'ATOM 2 C CA . THR A 1 2.450 2.000 3.000 95.40 1',
    'ATOM 3 C C . THR A 1 3.100 3.400 3.000 95.40 1',
    'ATOM 4 O O . THR A 1 4.300 3.550 3.000 95.40 1',
  ].join('\n');

  const { atoms, ssMap } = parseCifText(sampleCif);
  assert.equal(atoms.length, 4);
  const residues = buildResiduesFromAtoms(atoms, ssMap);
  assert.equal(residues.length, 1);
  assert.equal(residues[0].resName, 'THR');
  assert.equal(residues[0].oneLetter, 'T');
  const cov = computeCovalentBonds(atoms);
  assert.ok(cov.length >= 3);
  const analysis = analyzeAlphaFoldConfidence(residues, 'ALPHAFOLD_DB', 'AF-TEST-F1');
  assert.equal(analysis.fractionVeryHigh, 100);
});
