#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/phmodel.h>
#include <openbabel/elements.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>
#include <openbabel/bond.h>
#include <openbabel/generic.h>
#include <openbabel/forcefield.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
using namespace OpenBabel;

void test_Fix1912_PDBReading()
{
  // Reading from a PDB file should set the residues
  // and mark chains as perceived
  OBMolPtr mol = OBTestUtil::ReadFile("00T_ideal_het.pdb");
  OB_ASSERT(mol->HasChainsPerceived());
  OBAtom* atom = mol->GetAtom(1);
  OBResidue* res = atom->GetResidue();
  OB_REQUIRE(res != nullptr);
  OB_COMPARE(res->GetAtomID(atom), " N19");
  OB_COMPARE(res->GetChain(), 'A');
}

std::string remove_slashr(const char* smi)
{
  // Remove \r if present to normalise across platforms
  std::string ans;
  const char *p = smi;
  while (*p) {
    if (*p != '\r')
      ans += *p;
    p++;
  }
  return ans;
}

struct CdxData {
  const char* fname;
  const char* smi;
};

// Some basic reading of ChemDraw files
// Note that we don't correctly read radicals - TODO
// Also, converting ChemDraw doesn't work with the Read() interface, only Convert()
void test_ChemDraw_Basic()
{
  static const CdxData cdxData[] = {
    { "ethanol.cdx", "CCO\t\n" },
    // cyclohexane -> benzene reaction, plus another cyclohexane drawn on its own
    { "molrxnmix.cdx", "C1CCCCC1>>c1ccccc1\t\nC1CCCCC1\t\n" },
    { "MeCN.cdx", "CC#N\t\n"}
  };

  ios_base::openmode imode = ios_base::in | ios_base::binary;
  unsigned int size = sizeof(cdxData) / sizeof(CdxData);
  OBConversion conv;
  OB_REQUIRE(conv.SetInAndOutFormats("cdx", "smi"));
  std::stringstream outs;
  conv.SetOutStream(&outs);

  for (int i=0; i<size; ++i) {
    std::string fname = OBTestUtil::GetFilename(cdxData[i].fname);
    std::ifstream ifs(fname.c_str(), imode);
    OB_REQUIRE(ifs.good());
    conv.SetInStream(&ifs);
    outs.str("");
    conv.Convert();
    std::string out = outs.str();
    OB_COMPARE(remove_slashr(out.c_str()), cdxData[i].smi);
  }
}

// Very basic test for cdxml
void test_ChemDraw_XML_Basic()
{
  static const CdxData cdxmlData[] = {
      {"methanol.cdxml", "CO\t8\n"}};

  ios_base::openmode imode = ios_base::in;
  unsigned int size = sizeof(cdxmlData) / sizeof(CdxData);
  OBConversion conv;
  OB_REQUIRE(conv.SetInAndOutFormats("cdxml", "smi"));
  std::stringstream outs;
  conv.SetOutStream(&outs);

  for (int i = 0; i < size; ++i)
  {
    std::string fname = OBTestUtil::GetFilename(cdxmlData[i].fname);
    std::ifstream ifs(fname.c_str(), imode);
    OB_REQUIRE(ifs.good());
    conv.SetInStream(&ifs);
    outs.str("");
    conv.Convert();
    std::string out = outs.str();
    OB_COMPARE(remove_slashr(out.c_str()), cdxmlData[i].smi);
  }
}

// A basic test of functionality
void test_OBChemTsfm()
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("smi");
  
  // Notes to self: Problems with OBChemTsfm:
  // tsfm.Init("Cl[C:1]-[C:2]", "[C:1]=[C:2]"); // TODO: Need to change the API to take const char
  // Init should wipe the state so that OBChemTsfm can safely be reused

  conv.ReadString(&mol, "NCCBr");
  OBChemTsfm tsfm;
  std::string start("[N:1]-C-C");
  std::string end("[N+:1]");
  tsfm.Init(start, end);
  tsfm.Apply(mol);
  std::string out = conv.WriteString(&mol, true);
  OB_COMPARE(out, "[NH3+]CCBr");

  conv.ReadString(&mol, "ClCCBr");
  start = "Cl[C:1]-[C:2]";
  end = "[C:1]=[C:2]";
  OBChemTsfm b;
  b.Init(start, end);
  b.Apply(mol);
  out = conv.WriteString(&mol, true);
  OB_COMPARE(out, "ClC=CBr");

  conv.ReadString(&mol, "ClC(=O)[O]");
  start = "[#6]-[OD1:1]";
  end = "[#6]-[O-1:1]";
  OBChemTsfm c;
  c.Init(start, end);
  c.Apply(mol);
  out = conv.WriteString(&mol, true);
  OB_COMPARE(out, "ClC(=O)[O-]");

  conv.ReadString(&mol, "Cl[C]CBr");
  start = "Cl[C:1]-[C:2]";
  end = "[C:1]=[C:2]";
  OBChemTsfm d;
  d.Init(start, end);
  d.Apply(mol);
  out = conv.WriteString(&mol, true);
  OB_COMPARE(out, "Cl[C]=CBr");
}

// Open Babel was previously disappearing triple bonds when provided with SMILES
// containing a triple bond in an aromatic ring
void test_AromaticTripleBond()
{
  const char* smis[] = {
    "CCc1n[c]#[c]n1CC2CC(C(=O)O2)(c3ccccc3)c4ccccc4",
    "CCc1nc#cn1CC2CC(C(=O)O2)(c3ccccc3)c4ccccc4" };
  
  OBConversion conv;
  conv.SetInFormat("smi");

  for(int i=0; i<2; ++i) {
    OBMol mol;
    conv.ReadString(&mol, smis[i]);
    bool hasTripleBond = false;
    FOR_BONDS_OF_MOL(bond, mol) {
      if (bond->GetBondOrder()==3)
        hasTripleBond = true;
    }
    OB_ASSERT(hasTripleBond);
  }
}

// A segfault was occuring when a Universal SMILES was output after an InChIfied SMILES.
// This was due to short-circuit caching of InChIs on reading. The fix was to limit
// the situations when the cached value was used, but also to delete the cached value
// in this particular instance.
void test_Issue135_UniversalSmiles()
{
  // Test writing U smiles after I smiles
  OBConversion conv;
  conv.SetInFormat("smi");
  OBMol mol;
  conv.ReadString(&mol, "C(=O)([O-])C(=O)O");
  conv.SetOutFormat("smi");
  conv.SetOptions("I", OBConversion::OUTOPTIONS);
  std::string res = conv.WriteString(&mol, true);
  OB_COMPARE(res, "C(=O)(C(=O)O)[O-]");
  conv.SetOptions("U", OBConversion::OUTOPTIONS);
  res = conv.WriteString(&mol, true);
  OB_COMPARE(res, "C(=O)(C(=O)[O-])O");
}

// Reading an InChI and then adding hydrogens messed up the structure
void test_Issue134_InChI_addH()
{
  OBConversion conv;
  conv.SetInFormat("inchi");
  OBMol mol;
  conv.ReadString(&mol, "InChI=1S/C2H7NO/c1-2(3)4/h2,4H,3H2,1H3/t2-/m0/s1");
  OB_ASSERT(!mol.HasData(OBGenericDataType::VirtualBondData));
  mol.AddHydrogens();
  conv.SetOutFormat("smi");
  std::string res = conv.WriteString(&mol, true);
  OB_COMPARE(res, "C[C@@H](N)O");
}

// Delete hydrogens should not remove charged or isotopic hydrogens or [H][H] or [Cu][H][Cu]
// or hydrogens with assigned atom classes
void test_Issue178_DeleteHydrogens()
{
  OBConversion conv;
  conv.SetInFormat("smi");
  OBMol mol;
  // Test DeleteHydrogens() and DeleteNonPolarHydrogens()
  static const char *smi[] = { "C[H]", "[H][H]", "C[1H]", "C[H]C", "C[H+]" };
  int numHs[] = { 0, 2, 1, 1, 1 };
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 2; ++j) {
      conv.ReadString(&mol, smi[i]);
      if (j == 0)
        mol.DeleteHydrogens();
      else
        mol.DeleteNonPolarHydrogens();
      int myNumHs = 0;
      FOR_ATOMS_OF_MOL(atom, mol)
        if (atom->GetAtomicNum() == OBElements::Hydrogen)
          myNumHs++;
      OB_COMPARE(myNumHs, numHs[i]);
    }
  }
  // Test DeletePolarHydrogens()
  static const char *smiB[] = { "N[H]", "[H][H]", "N[1H]", "N[H]C", "N[H+]" };
  int numHsB[] = { 0, 2, 1, 1, 1 };
  for (int i = 0; i < 5; ++i) {
    conv.ReadString(&mol, smiB[i]);
    mol.DeletePolarHydrogens();
    int myNumHs = 0;
    FOR_ATOMS_OF_MOL(atom, mol)
      if (atom->GetAtomicNum() == OBElements::Hydrogen)
        myNumHs++;
    OB_COMPARE(myNumHs, numHsB[i]);
  }
  // Test atom class
  // Currently, the SMILES parser does not retain atom classes for hydrogens on reading so...
  conv.ReadString(&mol, "C[H]");
  OBPairInteger *ac = new OBPairInteger();
  ac->SetAttribute("Atom Class");
  ac->SetValue(99);
  mol.GetAtom(2)->SetData(ac); // Assign the hydrogen (atom 2) a class of 99
  mol.DeleteHydrogens();
  int myNumHs = 0;
  FOR_ATOMS_OF_MOL(atom, mol)
    if (atom->GetAtomicNum() == OBElements::Hydrogen)
      myNumHs++;
  OB_COMPARE(myNumHs, 1);
}

void test_Issue305_NumRotors()
{
  OBMolPtr mol = OBTestUtil::ReadFile("regressiontest_numrotors.mol");
  OB_COMPARE(mol->NumRotors(), 9); // was returning 4
}

void test_PR329_Molfile_RGroups()
{
  // There are several way to get an R Group into a mol file
  // 1. The existing use of atom maps on dummy atoms in SMILES
  OBConversion conv;
  OBMol mol;
  conv.SetInAndOutFormats("smi", "mol");
  conv.ReadString(&mol, "C([*:1]CO[*:2]");
  obErrorLog.SetOutputLevel(obError); // avoid warning about no 2D or 3D coords
  std::string molfileWithRGP = conv.WriteString(&mol);
  obErrorLog.SetOutputLevel(obWarning);
  OB_ASSERT( molfileWithRGP.find("R#") != std::string::npos );
  OB_ASSERT( molfileWithRGP.find("M  RGP  2   2   1   5   2") != std::string::npos); // i.e. atom 2 is labelled R1, atom 5 is labelled R2
  // Check negative case
  conv.ReadString(&mol, "C([*]CO[*]");
  std::string molfileb = conv.WriteString(&mol);
  OB_ASSERT( molfileb.find("R#") == std::string::npos );
  OB_ASSERT( molfileb.find("M  RGP") == std::string::npos);

  // 2. By reading a molfile that use the R#, RGP notation
  conv.SetInAndOutFormats("mol", "mol");
  conv.ReadString(&mol, molfileWithRGP);
  molfileb = conv.WriteString(&mol);
  OB_ASSERT( molfileb.find("R#") != std::string::npos );
  OB_ASSERT( molfileb.find("M  RGP  2   2   1   5   2") != std::string::npos); // i.e. atom 2 is labelled R1, atom 5 is labelled R2

  // 3. By reading a molfile that specifies the atom alias as Rn, where n is an integer
  std::string molfileWithAlias = "\n"
" OpenBabel07211621152D\n"
"\n"
"  2  1  0  0  0  0  0  0  0  0999 V2000\n"
"    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
"    0.0000    1.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n"
"  1  2  1  0  0  0  0\n"
"A    2\n"
"R1\n"
"M  END";
  conv.SetInAndOutFormats("mol", "mol");
  conv.ReadString(&mol, molfileWithAlias);
  std::string molfile = conv.WriteString(&mol);
  OB_ASSERT( molfile.find("R#") != std::string::npos );
  OB_ASSERT( molfile.find("M  RGP  1   2   1") != std::string::npos); // i.e. atom 2 is labelled R1
  // Check negative case
  molfileWithAlias = "\n"
" OpenBabel07211621152D\n"
"\n"
"  2  1  0  0  0  0  0  0  0  0999 V2000\n"
"    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
"    0.0000    1.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n"
"  1  2  1  0  0  0  0\n"
"A    2\n"
"R\n"
"M  END";
  conv.SetInAndOutFormats("mol", "mol");
  obErrorLog.SetOutputLevel(obError); // avoid warning "Alias R was not chemically interpreted"
  conv.ReadString(&mol, molfileWithAlias);
  obErrorLog.SetOutputLevel(obWarning);
  molfile = conv.WriteString(&mol);
  OB_ASSERT( molfile.find("R#") == std::string::npos );
  OB_ASSERT( molfile.find("M  RGP") == std::string::npos);

  // 4. By reading a molfile that specifies the element name as R1, etc.
  std::string molfileWithRGroupElementName = "\n"
" OpenBabel07211621152D\n"
"\n"
"  2  1  0  0  0  0  0  0  0  0999 V2000\n"
"    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
"    0.0000    1.0000    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0\n"
"  1  2  1  0  0  0  0\n"
"M  END";
  conv.SetInAndOutFormats("mol", "mol");
  conv.ReadString(&mol, molfileWithRGroupElementName);
  molfile = conv.WriteString(&mol);
  OB_ASSERT( molfile.find("R#") != std::string::npos );
  OB_ASSERT( molfile.find("M  RGP  1   2   1") != std::string::npos); // i.e. atom 2 is labelled R1
}

struct SmilesData {
  const char* inp;
  const char* out;
  const char* out_hoption; // if you specify -xh, i.e. "output explicit Hydrogens as such"
  const char* out_addh_hoption; // if you add hydrogens and then specify "-xh"
  const char* out_soption; // if you specify -xs, create SMARTS for substructure searching
};

void test_SMILES_Valence()
{
  static const SmilesData smilesData[] = {
    { "[H]", "", "", "", "[#1]" }, // SMARTS perhaps should be [*H], but not implemented at the moment
    { "[H][H]", "", "", "","[#1][#1]" },
    { "[HH]", "", "", "[H][H]", "[#1]" },
    { "C", "", "", "C([H])([H])([H])[H]", "" },
    { "[C]", "", "", "", "C" },
    { "[CH]", "", "", "[C][H]", "C" },
    { "[CH3]", "", "", "[C]([H])([H])[H]", "C" },
    { "[CH4]", "C", "C", "C([H])([H])([H])[H]", "C" },
    { "[CH5]", "", "", "[C]([H])([H])([H])([H])[H]", "C" },
    { "C[H]", "C", "", "C([H])([H])([H])[H]", "[C!H0]" },
    { "[C][H]", "[CH]", "", "", "[C!H0]" },
    { "[CH3][H]", "C", "C[H]", "C([H])([H])([H])[H]", "[C!H0]" },
    { "[CH2]([H])[H]", "C", "C([H])[H]", "C([H])([H])([H])[H]", "[C!H0!H1]" },
    { "[U][H]", "[UH]", "", "", "[U!H0]" },
    { "[UH2]", "", "", "[U]([H])[H]", "[U]" },
    { "[C@@H](Br)(Cl)I", "", "", "[C@](Br)(Cl)(I)[H]", "[C](Br)(Cl)I" }, // Note: if OB supported it, SMARTS should be [C@@?](Br)(Cl)I
    { "Br[C@@H](Cl)I", "", "", "Br[C@@](Cl)(I)[H]", "Br[C](Cl)I" }, // Note: if OB supported it, SMARTS should be Br[C@@?](Cl)I
    { "[C@@](F)(Br)(Cl)I", "", "", "", "" },
    { "F[C@@](Br)(Cl)I", "", "", "", "" },
    { "[H][C@@](Br)(Cl)I", "[C@@H](Br)(Cl)I", "", "", "[C@@H](Br)(Cl)I" },
    { "C[H:1]", "C", "C[H]", "C([H])([H])([H])[H]", "[C!H0]" }, // atom class only shown with -xa
    { "C[2H]", "", "", "C([2H])([H])([H])[H]", "C[2#1]" },
    { "c1ccccc1", "", "", "c1(c(c(c(c(c1[H])[H])[H])[H])[H])[H]", "" },
    { "c1cnccc1", "", "", "c1(c(nc(c(c1[H])[H])[H])[H])[H]", "" },
    { "c1c[nH]cc1", "", "", "c1(c(n(c(c1[H])[H])[H])[H])[H]", "c1cncc1" },
    { "F[I]F", "", "", "", "FIF" },
  };
  unsigned int size = (unsigned int)(sizeof(smilesData) / sizeof(smilesData[0]));
  for (unsigned int rep = 0; rep < 4; ++rep) {
    printf("Rep: %d\n", rep);
    OBConversion conv;
    OB_ASSERT(conv.SetInAndOutFormats("smi", "smi"));
    switch (rep) {
    case 1: case 2: conv.SetOptions("h", conv.OUTOPTIONS); break;
    case 3: conv.SetOptions("s", conv.OUTOPTIONS); break;
    }
    for (unsigned int i = 0; i < size; ++i) {
      OBMol mol;
      OB_ASSERT(conv.ReadString(&mol, smilesData[i].inp));
      if (rep == 2)
        mol.AddHydrogens();
      std::string out = conv.WriteString(&mol, true);
      const char* mout;
      switch (rep) {
      case 0: mout = smilesData[i].out; break;
      case 1: mout = smilesData[i].out_hoption; break;
      case 2: mout = smilesData[i].out_addh_hoption; break;
      case 3: mout = smilesData[i].out_soption; break;
      }
      std::string ans = mout[0] ? mout : smilesData[i].inp;
      printf("  %d %s --> %s (%s)\n", i, smilesData[i].inp, ans.c_str(), out.c_str());
      OB_COMPARE(out, ans);
    }
  }

  OBConversion conv;
  OB_ASSERT(conv.SetInAndOutFormats("smi", "smi"));
  conv.SetOptions("ah", conv.OUTOPTIONS); // write out alias explicitly
  OBMol mol;
  conv.ReadString(&mol, "C[H:1]");
  OB_COMPARE(conv.WriteString(&mol, true), "C[H:1]");
}

//make sure insertion code gets copied (it wasn't)
void test_insertioncode() 
{
  const char* pdb = "ATOM    266  HB2 ASP L  14      -2.604   8.021  19.867  1.00  0.00           H\n\
ATOM    267  CG  ASP L  14      -2.280   6.992  21.697  1.00 18.10           C\n\
ATOM    268  OD1 ASP L  14      -1.109   7.431  21.698  1.00 18.97           O\n\
ATOM    269  OD2 ASP L  14      -2.735   6.263  22.603  1.00 19.18           O\n\
ATOM    270  N   LYS L  14A     -5.804   6.060  21.469  1.00 20.85           N\n\
ATOM    271  H   LYS L  14A     -5.589   5.759  20.497  1.00  0.00           H\n\
ATOM    272  CA  LYS L  14A     -6.654   5.209  22.296  1.00 22.86           C\n\
ATOM    273  HA  LYS L  14A     -7.392   5.923  22.662  1.00  0.00           H\n\
ATOM    274  C   LYS L  14A     -6.108   4.607  23.591  1.00 21.70           C\n\
ATOM    275  O   LYS L  14A     -6.892   4.228  24.455  1.00 21.72           O\n";

    OBConversion conv;
    OB_ASSERT(conv.SetInAndOutFormats("pdb", "pdb"));
    OBMol mol;
    conv.ReadString(&mol, pdb);
    OBMol mol2;
    mol2 = mol;
    char i = mol2.GetResidue(1)->GetInsertionCode();
    OB_COMPARE(i, 'A');
}

//make sure icode is read by pdbqt
void test_insertioncode_pdbqt() 
{
  const char* pdb = "ATOM    266  HB2 ASP L  14      -2.604   8.021  19.867  1.00  0.00           H\n\
ATOM    267  CG  ASP L  14      -2.280   6.992  21.697  1.00 18.10           C\n\
ATOM    268  OD1 ASP L  14      -1.109   7.431  21.698  1.00 18.97           O\n\
ATOM    269  OD2 ASP L  14      -2.735   6.263  22.603  1.00 19.18           O\n\
ATOM    270  N   LYS L  14A     -5.804   6.060  21.469  1.00 20.85           N\n\
ATOM    271  H   LYS L  14A     -5.589   5.759  20.497  1.00  0.00           H\n\
ATOM    272  CA  LYS L  14A     -6.654   5.209  22.296  1.00 22.86           C\n\
ATOM    273  HA  LYS L  14A     -7.392   5.923  22.662  1.00  0.00           H\n\
ATOM    274  C   LYS L  14A     -6.108   4.607  23.591  1.00 21.70           C\n\
ATOM    275  O   LYS L  14A     -6.892   4.228  24.455  1.00 21.72           O\n";

    OBConversion conv;
    OB_ASSERT(conv.SetInAndOutFormats("pdbqt", "pdbqt"));
    OBMol mol;
    conv.ReadString(&mol, pdb);
    OBMol mol2;
    mol2 = mol;
    char i = mol2.GetResidue(1)->GetInsertionCode();
    OB_COMPARE(i, 'A');
}

// https://github.com/openbabel/openbabel/issues/1794
void test_github_issue_1794()
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "CC[2H]");

  OBForceField* pFF = OBForceField::FindForceField("UFF");
  OB_REQUIRE(pFF);

  OB_ASSERT(pFF->Setup(mol));
}

void test_github_issue_2111_impl(const std::string &smiles)
{
  OBMol mol;
  OBConversion conv;
  conv.SetInAndOutFormats("smi", "inchi");

  conv.ReadString(&mol, smiles);
  mol.DeleteHydrogens();
  conv.WriteString(&mol);
}

// https://github.com/openbabel/openbabel/issues/2111
void test_github_issue_2111()
{
  test_github_issue_2111_impl("[H][C@@H](I)F"); // tetrahedral with 2 implicit refs
  test_github_issue_2111_impl("[H]/N=C/F"); // cis/trans with 2 implicit refs on the left
  test_github_issue_2111_impl("F/N=C/[H]"); // cis/trans with 2 implicit refs on the right

  //
  // Tetrahedral
  //

  // example from bug report
  test_github_issue_2111_impl("[C@@H]12C[C@H](OCC3=CC=CC=C3)[C@@H](COCC3=CC=CC=C3)[C@]1([H])O2");

  // implicit ref in all positions
  test_github_issue_2111_impl("[H][C@](C)(F)I");
  test_github_issue_2111_impl("C[C@]([H])(F)I");
  test_github_issue_2111_impl("C[C@](F)([H])I");
  test_github_issue_2111_impl("C[C@](F)(I)[H]");

  //
  // Cis/Trans
  //

  // implicit ref in all positions
  test_github_issue_2111_impl("[H]/N(C)=C/F");
  test_github_issue_2111_impl("C/N([H])=C/F");
  test_github_issue_2111_impl("F/N=C(/[H])C");
  test_github_issue_2111_impl("F/N=C(/C)[H]");
}

void test_github_issue_2428_data_in_png()
{
  ios_base::openmode imode = ios_base::in | ios_base::binary;
  OBConversion conv;
  OB_REQUIRE(conv.SetInAndOutFormats("png", "can"));
  std::stringstream outs;
  conv.SetOutStream(&outs);

  std::string fname = OBTestUtil::GetFilename("pyridine.png");
  std::ifstream ifs(fname.c_str(), imode);
  OB_REQUIRE(ifs.good());
  conv.SetInStream(&ifs);
  outs.str("");
  conv.Convert();
  std::string out = outs.str();
  OB_COMPARE(remove_slashr(out.c_str()), "c1cccnc1\t\n");
}

void test_bad_bondorders() //make sure N in aromatic rings get right bond orders
{
    //from PDB 1U71
  const char* pdb = "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE MXA A 187\n\
HET    MXA  A 187      24                                                       \n\
HETNAM     MXA 6-(2,5-DIMETHOXY-BENZYL)-5-METHYL-PYRIDO[2,3-                    \n\
HETNAM   2 MXA  D]PYRIMIDINE-2,4-DIAMINE                                        \n\
FORMUL   4  MXA    C17 H19 N5 O2                                                \n\
HETATM 1517  N1  MXA A 187      28.796  13.139  -4.319  0.50  2.21           N  \n\
HETATM 1518  C2  MXA A 187      27.596  12.644  -4.455  0.50  2.00           C  \n\
HETATM 1519  N2  MXA A 187      27.196  12.327  -5.700  0.50  2.00           N  \n\
HETATM 1520  N3  MXA A 187      26.752  12.365  -3.424  0.50  2.00           N  \n\
HETATM 1521  C4  MXA A 187      27.102  12.613  -2.212  0.50  2.38           C  \n\
HETATM 1522  N4  MXA A 187      26.137  12.302  -1.332  0.50  2.06           N  \n\
HETATM 1523  C4A MXA A 187      28.375  13.227  -1.915  0.50  3.60           C  \n\
HETATM 1524  C5  MXA A 187      28.903  13.607  -0.674  0.50  5.54           C  \n\
HETATM 1525  C5M MXA A 187      28.179  13.440   0.652  0.50  5.83           C  \n\
HETATM 1526  C6  MXA A 187      30.171  14.172  -0.645  0.50  7.15           C  \n\
HETATM 1527  C7  MXA A 187      30.848  14.357  -1.861  0.50  5.90           C  \n\
HETATM 1528  N8  MXA A 187      30.407  14.021  -3.042  0.50  3.72           N  \n\
HETATM 1529  C8A MXA A 187      29.188  13.455  -3.076  0.50  3.31           C  \n\
HETATM 1530  C9  MXA A 187      30.823  14.676   0.605  0.50 10.20           C  \n\
HETATM 1531  C1' MXA A 187      32.247  14.320   0.954  0.50 12.53           C  \n\
HETATM 1532  C2' MXA A 187      33.277  15.250   0.897  0.50 13.42           C  \n\
HETATM 1533  O2' MXA A 187      32.883  16.513   0.499  0.50 13.85           O  \n\
HETATM 1534  C21 MXA A 187      33.905  17.478   0.251  0.50 13.61           C  \n\
HETATM 1535  C3' MXA A 187      34.547  14.874   1.267  0.50 14.29           C  \n\
HETATM 1536  C4' MXA A 187      34.808  13.570   1.706  0.50 14.41           C  \n\
HETATM 1537  C5' MXA A 187      33.819  12.674   1.740  0.50 13.91           C  \n\
HETATM 1538  O5' MXA A 187      33.800  11.326   2.099  0.50 16.10           O  \n\
HETATM 1539  C51 MXA A 187      32.462  10.825   2.070  0.50 15.86           C  \n\
HETATM 1540  C6' MXA A 187      32.533  13.060   1.396  0.50 13.52           C  ";

    OBConversion conv;
    OB_ASSERT(conv.SetInAndOutFormats("pdb", "sdf"));
    OBMol mol;
    conv.ReadString(&mol, pdb);
    //no nitrogens should have 4 bonds
    for (OBAtomIterator it=mol.BeginAtoms(); it != mol.EndAtoms(); it++)
    {
        OBAtom *a = *it;
        if(a->GetAtomicNum() == 7) {
        	OB_ASSERT((a->GetTotalValence() < 4));
        }
    }
}
    
void test_SegCopySubstructure()
{
  // Invalid memory access (atom->GetIdx()) detected in valgrind and sometimes
  // triggering a sefault.
  OBConversion conv;
  OB_ASSERT(conv.SetInFormat("smi"));
  OBMol mol;
  std::string smi = "C[C@@H]1CO1";
  OB_ASSERT(conv.ReadString(&mol, smi));

  OBBitVec atomsToCopy;
  atomsToCopy.Clear();
  atomsToCopy.SetBitOn(2);
  atomsToCopy.SetBitOn(3);
  atomsToCopy.SetBitOn(4);
  OBBitVec bondsToExclude;
  bondsToExclude.Clear();
  bondsToExclude.SetBitOn(0);

  OBMol copy;
  std::vector<unsigned int> atomorder;
  atomorder.clear();
  bool ok = mol.CopySubstructure(copy, &atomsToCopy, &bondsToExclude, 2, &atomorder);
  OB_ASSERT(ok);
  OB_COMPARE(4, copy.NumAtoms());
  OB_COMPARE(4, copy.NumBonds());
}

int regressionstest(int argc, char* argv[])
{
  int defaultchoice = 1;
  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }
  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif  

  switch(choice) {
  case 1:
    test_Issue135_UniversalSmiles();
    break;
  case 2:
    test_SegCopySubstructure();
    break;
  case 221:
    test_Issue134_InChI_addH();
    break;
  case 222:
    test_Issue178_DeleteHydrogens();
    break;
  case 223:
    test_Issue305_NumRotors();
    break;
  case 224:
    test_PR329_Molfile_RGroups();
    break;
  case 225:
    test_AromaticTripleBond();
    break;
  case 226:
    test_SMILES_Valence();
    break;
  case 227:
    test_OBChemTsfm();
    break;
  case 228:
    test_ChemDraw_Basic();
    break;
  case 229:
    test_ChemDraw_XML_Basic();
    break;
  case 240:
    test_Fix1912_PDBReading();
    break;
  case 241:
    test_insertioncode();
    break;
  case 242:
    test_insertioncode_pdbqt();
    break;
  case 1794:
    test_github_issue_1794();
    break;
  case 2111:
    test_github_issue_2111();
    break;
  case 2428:
    test_github_issue_2428_data_in_png();
    break;
  case 2384:
	  test_bad_bondorders();
	  break;
  //case N:
  //  YOUR_TEST_HERE();
  //  Remember to update CMakeLists.txt with the number of your test
  //  break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}

