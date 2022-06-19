/**********************************************************************
Copyright (C) 2004 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/


#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/oberror.h>
#include <openbabel/descriptor.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/residue.h>
#include <ctime>
#include <vector>
#include <iomanip>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif

using namespace std;

namespace OpenBabel
{

    class CTLDFormat : public OBMoleculeFormat
        // Derive directly from OBFormat for objects which are not molecules.
    {
    public:

        string rdfType;
        string writeCompound;
        string writeHasAtom;
        string writeAtom;
        string writeatom;
        string writexCoordinate;
        string writeyCoordinate;
        string writezCoordinate;
        string writeXMLdecimal;
        string writeXMLstring;
        string writeEndPoint;
        string writeHasBond;
        string writeBond;
        string writeFirstAtom;
        string writeSecondAtom;
        string writeType;
        string writeSingle;
        string writeDouble;
        string writeStereo;
        string writeline;
        string writeColon;
        string writePrefix;
        string enrichFormula;
        string enrichMolecularWeight;
        string enrichExactMass;
        string enrichLogP;
        string enrichTpsa;
        string enrichMolecularRefractivity;
        string enrichSmiles;
        string enrichInchi;
        string enrichInchiKey;
        //Register this format type ID in the constructor


        /* The first line of the description should be a brief identifier, <40 chars, because
         it is used in dropdown lists, etc. in some user interfaces. The rest is optional.

           Describe any format specific options here. This text is parsed to provide
           checkboxes, etc for the GUI (for details click the control menu),
           so please try to keep to a similar form.

           Write options are the most common, and the "Write" is optional.
           The option f takes a text parameter, so that it is essential that the option
           is registered in the constructor of the class.
           Finish the options with a blank line as shown, if there are more than one
           group of options, or if there are further comments after them.
        */
        virtual const char* Description() //required
        {
            return
                "CTLD format\n"
                "CT-LD distinguishes three main fragments that contain triples of RDF.\n\n"
                "The compound triples contain information about the molecule, e.g., molar refractivity (MR), polar surface area (PSA), logarithm of the octanol–water partition coefficient (LogP) etc. Here, it is possible to embed line notation formats, e.g., SMILES, InChI, InChIKey, etc. \n\n"
                "Atom triples define the atomic symbol and any mass difference, charge, stereochemistry, and associated hydrogens for each atom. \n\n"
                "Bond triples define the two atoms connected by the bond, the bond type, and any bond stereochemistry and topology (chain or ring properties) for each bond.\n\n"
                "Ontology http://ii.uwb.edu.pl/ctld \n\n"
                "Write Options e.g. -xE \n"
                "	E enrich informations\n"	

                "Read Options e.g. -aS\n"
                "	S  Sorted file\n"
                ;
        };

        //Optional URL where the file format is specified
        virtual const char* SpecificationURL() { return "http://ii.uwb.edu.pl/ctld"; }

        //Optional
        virtual const char* GetMIMEType()
        {
            return "chemical/x-ctld";
        };


        /* Flags() can return be any of the following combined by |
             or be omitted if none apply
           NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
           READBINARY  WRITEBINARY  READXML  ZEROATOMSOK*/
           //   virtual unsigned int Flags()
           //   {
           //       return READONEONLY;
           //   };

               /* This optional function is for formats which can contain more than one
                  molecule. It is used to quickly position the input stream after the nth
                  molecule without have to convert and discard all the n molecules.
                  See obconversion.cpp for details and mdlformat.cpp for an example.*/
                  // virtual int SkipObjects(int n, OBConversion* pConv)
                  // {
                  // 	return 0;
                  // };

                  ////////////////////////////////////////////////////
                /// Declarations for the "API" interface functions. Definitions are below
        virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
        virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    private:
        enum Parity {
            NotStereo, Clockwise, AntiClockwise, Unknown
        };
        void GetParity(OBMol& mol, map<OBAtom*, Parity>& parity);
        void EscapeRDF(string& str);
        int writeMolIterator = 0;
        bool isFirstRead = true;
        bool firstMolPopulation = true;
        int whatCompound = 0;
        bool endMolPopulation = false;
        bool prefixApplied = false;

        //Vectors to hold data
        std::vector<string> compoundWithIndex;


        std::vector<vector<string>> compoundAtomsVect;
        std::vector<vector<string>> compoundAtomsIndxVect;

        std::vector<vector<string>> AtomsInfoX;
        std::vector<vector<string>> AtomsInfoY;
        std::vector<vector<string>> AtomsInfoZ;
        std::vector<vector<string>> AtomsInfoAtom;


        std::vector<vector<string>> compoundBondsVect;

        std::vector<vector<string>> BondsInfoFirst;
        std::vector<vector<string>> BondsInfoSecond;
        std::vector<vector<string>> BondsInfoType;



        //Special delimiters
        string dotDelimiter = ".";
        string semicolonDelimiter = ";";

        //Const file variables needed to read
        string isCompound = ":Compound";
        string isChiral = ":chiral";
        string isAtom = ":Atom";
        string x_coordinate = ":xCoordinate";
        string y_coordinate = ":yCoordinate";
        string z_coordinate = ":zCoordinate";
        string atomName = " :atom";
        string isStereo = ":hasStereoConf";
        string delimiterNorm = "_";
        string delimiterLinkOpen = "<";
        string delimiterLinkClose = ">";
        string bondType = ":hasBondType";
        string isFirst = ":firstAtom";
        string isSecond = ":secondAtom";
        string isBond = ":Bond";
        string hasAtom = ":hasAtom";
        string hasBond = ":hasBond";

        string bondsType = ":type";
        string bondsTypeLn = "#type>";

        // For link type file
        string isCompoundLn = "#Compound>";
        string isChiralLn = "#chiral>";
        string isAtomLn = "#Atom>";
        string x_coordinateLn = "#xCoordinate>";
        string y_coordinateLn = "#yCoordinate>";
        string z_coordinateLn = "#zCoordinate>";
        string atomNameLn = "#atom>";
        string isStereoLn = "#hasStereoConf>";
        string bondTypeLn = "#hasBondType>";
        string isFirstLn = "#firstAtom>";
        string isSecondLn = "#secondAtom>";
        string isBondLn = "#Bond>";
        string hasAtomLn = "#hasAtom>";
        string hasBondLn = "#hasBond>";

        //ns#type
        string nsType = "ns#type>";


        string actualCompound;



        int atoms_count = 0; //Global atoms count
        int actualAtom = 0; //Current atom input
        int bonds_count = 0; // Global bonds count
        int actualBond = 0; // Current bond input;

    //Counters and variables for ns#type>
        int whichCompound = 0;
        string currentCompound;
        string currentAtom;
        string currentBond;

        int compIter = 0;

        int whatComp = 1;
        bool found = false;
        bool foundX = false;
        bool foundY = false;
        bool foundZ = false;
        bool foundAtm = false;
        bool foundBndFir = false;
        bool foundBndSec = false;
        bool foundBndType = false;

        int moleculeIterator = 0;
        int moleculeIteratorRead = 0;
        int atmNmIndx = 1;
        bool firstComp = false;
        // 	/* Add declarations for any local function or member variables used.
        // 	   Generally only a single instance of a format class is used. Keep this in
        // 	   mind if you employ member variables. */
    };
    ////////////////////////////////////////////////////

//Make an instance of the format class

    class TTLFormat : public CTLDFormat
    {
    public:
        //Register this format type ID
        TTLFormat()
        {
            OBConversion::RegisterFormat("ttl", this, "chemical/x-ctld-ttlfile");
            OBConversion::RegisterOptionParam("S", this);
            OBConversion::RegisterOptionParam("E", this, OBConversion::OUTOPTIONS);
        }

        virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv)
        {
            rdfType = "a";
            writeCompound = ":Compound";
            writeHasAtom = ":hasAtom";
            writeAtom = ":Atom";
            writeatom = ":atom";
            writexCoordinate = ":xCoordinate";
            writeyCoordinate = ":yCoordinate";
            writezCoordinate = ":zCoordinate";
            writeXMLdecimal = "";
            writeXMLstring = "";
            writeEndPoint = " .";
            writeHasBond = ":hasBond";
            writeBond = ":Bond";
            writeFirstAtom = ":firstAtom";
            writeSecondAtom = ":secondAtom";
            writeType = ":type";
            writeSingle = ":single";
            writeDouble = ":double";
            writeline = "\n";
            writeStereo = ":stereo";
            writeColon = "";
            writePrefix = "@prefix : <https://ii.uwb.edu.pl/ctld#> .";
            enrichFormula = ":formula";
            enrichMolecularWeight = ":molecularWeight";
            enrichExactMass = ":exactMass";
            enrichLogP = ":logP";
            enrichTpsa = ":tpsa";
            enrichMolecularRefractivity = ":molecularRefractivity";
            enrichSmiles = ":smiles";
            enrichInchi = ":inchi";
            enrichInchiKey = ":inchikey";
            
            return CTLDFormat::WriteMolecule(pOb, pConv);
        }
    };

    //Make an instance of the format class
    TTLFormat theTTLFormat;

    //*************************************
    class NTFormat : public CTLDFormat
    {
    public:
        NTFormat()
        {
            OBConversion::RegisterFormat("nt", this, "chemical/x-ctld-ntfile");
            OBConversion::RegisterOptionParam("S", this);
            OBConversion::RegisterOptionParam("E", this, OBConversion::OUTOPTIONS);
        }

        virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv)
        {
            rdfType = "<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>";
            writeCompound = "<https://ii.uwb.edu.pl/ctld#Compound>";
            writeHasAtom = "<https://ii.uwb.edu.pl/ctld#hasAtom>";
            writeAtom = "<https://ii.uwb.edu.pl/ctld#Atom>";
            writeatom = "<https://ii.uwb.edu.pl/ctld#atom>";
            writexCoordinate = "<https://ii.uwb.edu.pl/ctld#xCoordinate>";
            writeyCoordinate = "<https://ii.uwb.edu.pl/ctld#yCoordinate>";
            writezCoordinate = "<https://ii.uwb.edu.pl/ctld#zCoordinate>";
            writeXMLdecimal = "^^<https://www.w3.org/2001/XMLSchema#decimal>";
            writeXMLstring = "^^<https://www.w3.org/2001/XMLSchema#string>";
            writeEndPoint = " .";
            writeHasBond = "<https://ii.uwb.edu.pl/ctld#hasBond>";
            writeBond = "<https://ii.uwb.edu.pl/ctld#Bond>";
            writeFirstAtom = "<https://ii.uwb.edu.pl/ctld#firstAtom>";
            writeSecondAtom = "<https://ii.uwb.edu.pl/ctld#secondAtom>";
            writeType = "<https://ii.uwb.edu.pl/ctld#type>";
            writeSingle = "<https://ii.uwb.edu.pl/ctld#single>";
            writeDouble = "<https://ii.uwb.edu.pl/ctld#double>";
            writeStereo = "<https://ii.uwb.edu.pl/ctld#stereo>";
            writeline = "";
            writeColon = "\"";
            writePrefix = "";
            enrichFormula = "<https://ii.uwb.edu.pl/ctld#formula>";
            enrichMolecularWeight = "<https://ii.uwb.edu.pl/ctld#molecularWeight>";
            enrichExactMass = "<https://ii.uwb.edu.pl/ctld#exactMass>";
            enrichLogP = "<https://ii.uwb.edu.pl/ctld#logP>";
            enrichTpsa = "<https://ii.uwb.edu.pl/ctld#tpsa>";
            enrichMolecularRefractivity = "<https://ii.uwb.edu.pl/ctld#molecularRefractivity>";
            enrichSmiles = "<https://ii.uwb.edu.pl/ctld#smiles>";
            enrichInchi = "<https://ii.uwb.edu.pl/ctld#inchi>";
            enrichInchiKey = "<https://ii.uwb.edu.pl/ctld#inchikey>";

            return CTLDFormat::WriteMolecule(pOb, pConv);
        }
    };

    //Make an instance of the format class
    NTFormat theNTFormat;




    /////////////////////////////////////////////////////////////////

    bool CTLDFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
    {
        OBMol* pmol = pOb->CastAndClear<OBMol>();
        if (pmol == nullptr)
            return false;
        OBMol& mol = *pmol;
        stringstream errorMsg;
        istream& ifs = *pConv->GetInStream();

        if (!ifs.good() || ifs.peek() == EOF)
            return false;

        //Variables needed to preprocess the file.

        string line; //Line from file




        //If file was sucessfully opened 

        if (isFirstRead) {
            cerr << "Started processing file" << "\n";

            if (std::getline(ifs, line))
            {
                while (std::getline(ifs, line)) {

                    if ((line[0] == '_' && line.find(isCompound) != string::npos && line[line.length() - 1] == ';') || (line[0] == '_' && line.find(isCompoundLn) != string::npos && line[line.length() - 1] == ';') || (line[0] == '_' && line.find(isCompound) != string::npos && line[line.length() - 1] == '.') || (line[0] == '_' && line.find(isCompoundLn) != string::npos && line[line.length() - 1] == '.')) {

                        size_t x = line.find(" ", 0);

                        currentCompound = line.substr(0, x);



                        char last = line[line.length() - 1];
                        if (find(compoundWithIndex.begin(), compoundWithIndex.end(), currentCompound) == compoundWithIndex.end())
                        {
                            compoundWithIndex.push_back(currentCompound);
                            compoundAtomsVect.push_back(vector<string>());
                            compoundAtomsIndxVect.push_back(vector<string>());
                            compoundBondsVect.push_back(vector<string>());
                            AtomsInfoX.push_back(vector<string>());
                            AtomsInfoY.push_back(vector<string>());
                            AtomsInfoZ.push_back(vector<string>());
                            AtomsInfoAtom.push_back(vector<string>());
                            BondsInfoFirst.push_back(vector<string>());
                            BondsInfoSecond.push_back(vector<string>());
                            BondsInfoType.push_back(vector<string>());
                        }


                        if (last == '.') {
                            currentCompound = "";
                        }

                    }
                    if ((currentCompound.length() > 2 && line.find(hasAtom) != string::npos) || (currentCompound.length() > 2 && line.find(hasAtomLn) != string::npos)) {

                        size_t y = 0;
                        if (line.find(hasAtomLn) != string::npos) {
                            y = line.find(hasAtomLn);
                        }
                        else if (line.find(hasAtom) != string::npos) {
                            y = line.find(hasAtom);
                        }

                        char last = line[line.length() - 1];

                        for (int i = y; i < line.length(); i++) {
                            if (line[i] == '_') {

                                atoms_count++;

                                size_t v = line.find(" ", i);

                                string atmNm = line.substr(i, (v - i));


                                if (find(compoundWithIndex.begin(), compoundWithIndex.end(), currentCompound) == compoundWithIndex.end())
                                {
                                    compoundWithIndex.push_back(currentCompound);
                                    compoundAtomsVect.push_back(vector<string>());
                                    compoundAtomsIndxVect.push_back(vector<string>());
                                    compoundBondsVect.push_back(vector<string>());
                                    AtomsInfoX.push_back(vector<string>());
                                    AtomsInfoY.push_back(vector<string>());
                                    AtomsInfoZ.push_back(vector<string>());
                                    AtomsInfoAtom.push_back(vector<string>());
                                    BondsInfoFirst.push_back(vector<string>());
                                    BondsInfoSecond.push_back(vector<string>());
                                    BondsInfoType.push_back(vector<string>());
                                    whatComp = compoundWithIndex.size() - 1;
                                    compoundAtomsVect[whatComp].push_back(atmNm);
                                    compoundAtomsIndxVect[whatComp].push_back(atmNm);
                                    compoundAtomsIndxVect[whatComp].push_back(std::to_string(atmNmIndx));
                                    atmNmIndx++;

                                }
                                else {
                                    auto it = find(compoundWithIndex.begin(), compoundWithIndex.end(), currentCompound);
                                    whatComp = it - compoundWithIndex.begin();
                                    compoundAtomsVect[whatComp].push_back(atmNm);
                                    compoundAtomsIndxVect[whatComp].push_back(atmNm);
                                    compoundAtomsIndxVect[whatComp].push_back(std::to_string(atmNmIndx));
                                    atmNmIndx++;
                                }

                            }
                        }
                        if (last == '.') {
                            currentCompound = "";
                        }
                    }
                    if ((currentCompound.length() == 0 && line.find(hasAtom) != string::npos && line[0] == '_') || (currentCompound.length() == 0 && line.find(hasAtomLn) != string::npos && line[0] == '_')) {

                        size_t x = line.find(" ", 0);

                        currentCompound = line.substr(0, x);

                        size_t y = 0;
                        if (line.find(hasAtomLn) != string::npos) {
                            y = line.find(hasAtomLn);
                        }
                        else if (line.find(hasAtom) != string::npos) {
                            y = line.find(hasAtom);
                        }

                        char last = line[line.length() - 1];

                        for (int i = y; i < line.length(); i++) {
                            if (line[i] == '_') {

                                atoms_count++;

                                size_t v = line.find(" ", i);

                                string atmNm = line.substr(i, (v - i));


                                if (find(compoundWithIndex.begin(), compoundWithIndex.end(), currentCompound) == compoundWithIndex.end())
                                {
                                    compoundWithIndex.push_back(currentCompound);
                                    compoundAtomsVect.push_back(vector<string>());
                                    compoundAtomsIndxVect.push_back(vector<string>());
                                    compoundBondsVect.push_back(vector<string>());
                                    AtomsInfoX.push_back(vector<string>());
                                    AtomsInfoY.push_back(vector<string>());
                                    AtomsInfoZ.push_back(vector<string>());
                                    AtomsInfoAtom.push_back(vector<string>());
                                    BondsInfoFirst.push_back(vector<string>());
                                    BondsInfoSecond.push_back(vector<string>());
                                    BondsInfoType.push_back(vector<string>());
                                    whatComp = compoundWithIndex.size() - 1;
                                    compoundAtomsVect[whatComp].push_back(atmNm);
                                    compoundAtomsIndxVect[whatComp].push_back(atmNm);
                                    compoundAtomsIndxVect[whatComp].push_back(std::to_string(atmNmIndx));
                                    atmNmIndx++;

                                }
                                else {
                                    auto it = find(compoundWithIndex.begin(), compoundWithIndex.end(), currentCompound);
                                    whatComp = it - compoundWithIndex.begin();
                                    compoundAtomsVect[whatComp].push_back(atmNm);
                                    compoundAtomsIndxVect[whatComp].push_back(atmNm);
                                    compoundAtomsIndxVect[whatComp].push_back(std::to_string(atmNmIndx));
                                    atmNmIndx++;
                                }


                            }
                        }
                        if (last == '.') {
                            currentCompound = "";
                        }
                    }
                    if ((currentCompound.length() > 2 && line.find(hasBond) != string::npos) || (currentCompound.length() > 2 && line.find(hasBondLn) != string::npos)) {
                        \

                            size_t x = 0;
                        if (line.find(hasBondLn) != string::npos) {
                            x = line.find(hasBondLn);
                        }
                        else if (line.find(hasBond) != string::npos) {
                            x = line.find(hasBond);
                        }

                        char last = line[line.length() - 1];

                        for (int i = x; i < line.length(); i++) {
                            if (line[i] == '_') {
                                bonds_count++;


                                size_t y = 0;
                                if (line.find(hasBondLn) != string::npos) {
                                    y = line.find(hasBondLn);
                                }
                                else if (line.find(hasBond) != string::npos) {
                                    y = line.find(hasBond);
                                }
                                size_t c = line.find(" ", y + 1);
                                size_t v = line.find(" ", c + 1);
                                string bndNm = line.substr(c + 1, (v - c) - 1);


                                if (find(compoundWithIndex.begin(), compoundWithIndex.end(), currentCompound) == compoundWithIndex.end())
                                {
                                    compoundWithIndex.push_back(currentCompound);
                                    compoundAtomsVect.push_back(vector<string>());
                                    compoundAtomsIndxVect.push_back(vector<string>());
                                    compoundBondsVect.push_back(vector<string>());
                                    AtomsInfoX.push_back(vector<string>());
                                    AtomsInfoY.push_back(vector<string>());
                                    AtomsInfoZ.push_back(vector<string>());
                                    AtomsInfoAtom.push_back(vector<string>());
                                    BondsInfoFirst.push_back(vector<string>());
                                    BondsInfoSecond.push_back(vector<string>());
                                    BondsInfoType.push_back(vector<string>());
                                    whatComp = compoundWithIndex.size() - 1;
                                    compoundBondsVect[whatComp].push_back(bndNm);

                                }
                                else {
                                    auto it = find(compoundWithIndex.begin(), compoundWithIndex.end(), currentCompound);
                                    whatComp = it - compoundWithIndex.begin();
                                    compoundBondsVect[whatComp].push_back(bndNm);
                                }
                            }
                        }
                        if (last == '.') {
                            currentCompound = "";
                        }
                    }
                    if ((currentCompound.length() == 0 && line.find(hasBond) != string::npos) || (currentCompound.length() == 0 && line.find(hasBondLn) != string::npos)) {

                        size_t x = line.find(" ", 0);

                        currentCompound = line.substr(0, x);


                        size_t y = 0;
                        if (line.find(hasBondLn) != string::npos) {
                            y = line.find(hasBondLn);
                        }
                        else if (line.find(hasBond) != string::npos) {
                            y = line.find(hasBond);
                        }

                        char last = line[line.length() - 1];

                        for (int i = y; i < line.length(); i++) {
                            if (line[i] == '_') {
                                bonds_count++;

                                x = y;
                                size_t c;
                                if (line.find(hasBondLn) != string::npos) {
                                    c = line.find_first_not_of(" \t", x + hasBondLn.length());
                                }
                                else if (line.find(hasBond) != string::npos) {
                                    c = line.find_first_not_of(" \t", x + hasBond.length());
                                }


                                size_t v = line.find(" ", c + 1);
                                string bndNm = line.substr(c, (v - c));


                                if (find(compoundWithIndex.begin(), compoundWithIndex.end(), currentCompound) == compoundWithIndex.end())
                                {
                                    compoundWithIndex.push_back(currentCompound);
                                    compoundAtomsVect.push_back(vector<string>());
                                    compoundAtomsIndxVect.push_back(vector<string>());
                                    compoundBondsVect.push_back(vector<string>());
                                    AtomsInfoX.push_back(vector<string>());
                                    AtomsInfoY.push_back(vector<string>());
                                    AtomsInfoZ.push_back(vector<string>());
                                    AtomsInfoAtom.push_back(vector<string>());
                                    BondsInfoFirst.push_back(vector<string>());
                                    BondsInfoSecond.push_back(vector<string>());
                                    BondsInfoType.push_back(vector<string>());
                                    whatComp = compoundWithIndex.size() - 1;
                                    compoundBondsVect[whatComp].push_back(bndNm);

                                }
                                else {
                                    auto it = find(compoundWithIndex.begin(), compoundWithIndex.end(), currentCompound);
                                    whatComp = it - compoundWithIndex.begin();
                                    compoundBondsVect[whatComp].push_back(bndNm);
                                }
                            }
                        }
                        if (last == '.') {
                            currentCompound = "";
                        }
                    }

                }

                ifs.clear();
                ifs.seekg(0);


                while (std::getline(ifs, line)) {

                    //Reading all atoms and bonds info.
                    if (pConv->IsOption("S", OBConversion::INOPTIONS))
                    {
                        if ((line[0] == '_' && line.find(isCompound) != string::npos && line[line.length() - 1] == ';') || (line[0] == '_' && line.find(isCompoundLn) != string::npos && line[line.length() - 1] == ';') || (line[0] == '_' && line.find(isCompound) != string::npos && line[line.length() - 1] == '.') || (line[0] == '_' && line.find(isCompoundLn) != string::npos && line[line.length() - 1] == '.')) {

                            if (firstComp) {
                                moleculeIteratorRead++;
                            }
                            else {
                                firstComp = true;
                            }

                        }
                    }
                    //ATOMS INFOS
                    if ((line[0] == '_' && line.find(isAtom) != string::npos && line[line.length() - 1] == ';') || (line[0] == '_' && line.find(isAtomLn) != string::npos && line[line.length() - 1] == ';') || (line[0] == '_' && line.find(isAtom) != string::npos && line[line.length() - 1] == '.') || (line[0] == '_' && line.find(isAtomLn) != string::npos && line[line.length() - 1] == '.')) {

                        size_t x = line.find(" ", 0);
                        char last = line[line.length() - 1];

                        currentAtom = line.substr(0, x);


                        if (last == '.') {
                            currentAtom = "";
                        }

                    }
                    if ((line.find(x_coordinate) != string::npos && currentAtom.length() > 2) || (line.find(x_coordinateLn) != string::npos && currentAtom.length() > 2)) {

                        size_t x = 0;
                        if (line.find(x_coordinateLn) != string::npos) {
                            x = line.find(x_coordinateLn);
                        }
                        else if (line.find(x_coordinate) != string::npos) {
                            x = line.find(x_coordinate);
                        }

                        char last = line[line.length() - 1];

                        string coords;


                        if (line.find("\"") != string::npos && line.find(x_coordinateLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);


                        }
                        else if (line.find(x_coordinateLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + x_coordinateLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(x_coordinate) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(x_coordinate) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + x_coordinate.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {

                            for (int j = 0; j < compoundAtomsVect[moleculeIteratorRead].size(); j++) {
                                if (compoundAtomsVect[moleculeIteratorRead][j] == currentAtom) {
                                    found = true;
                                    AtomsInfoX[moleculeIteratorRead].push_back(currentAtom);
                                    AtomsInfoX[moleculeIteratorRead].push_back(coords);
                                    break;
                                }
                            }
                        }
                        else {

                            for (int i = 0; i < compoundAtomsVect.size(); i++) {
                                for (int j = 0; j < compoundAtomsVect[i].size(); j++) {

                                    if (compoundAtomsVect[i][j] == currentAtom) {
                                        found = true;
                                        AtomsInfoX[i].push_back(currentAtom);
                                        AtomsInfoX[i].push_back(coords);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;

                        if (last == '.') {
                            currentAtom = "";
                        }
                    }
                    if ((currentAtom.length() > 2 && line.find(y_coordinate) != string::npos) || (currentAtom.length() > 2 && line.find(y_coordinateLn) != string::npos)) {

                        size_t x = 0;
                        if (line.find(y_coordinateLn) != string::npos) {
                            x = line.find(y_coordinateLn);
                        }
                        else if (line.find(y_coordinate) != string::npos) {
                            x = line.find(y_coordinate);
                        }

                        char last = line[line.length() - 1];

                        string coords;

                        if (line.find("\"") != string::npos && line.find(y_coordinateLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(y_coordinateLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + y_coordinateLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(y_coordinate) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(y_coordinate) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + y_coordinate.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }


                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundAtomsVect[moleculeIteratorRead].size(); j++) {
                                if (compoundAtomsVect[moleculeIteratorRead][j] == currentAtom) {
                                    found = true;
                                    AtomsInfoY[moleculeIteratorRead].push_back(currentAtom);
                                    AtomsInfoY[moleculeIteratorRead].push_back(coords);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundAtomsVect.size(); i++) {
                                for (int j = 0; j < compoundAtomsVect[i].size(); j++) {

                                    if (compoundAtomsVect[i][j] == currentAtom) {
                                        found = true;
                                        AtomsInfoY[i].push_back(currentAtom);
                                        AtomsInfoY[i].push_back(coords);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;


                        if (last == '.') {
                            currentAtom = "";
                        }
                    }
                    if ((currentAtom.length() > 2 && line.find(z_coordinate) != string::npos) || (currentAtom.length() > 2 && line.find(z_coordinateLn) != string::npos)) {

                        size_t x = 0;
                        if (line.find(z_coordinateLn) != string::npos) {
                            x = line.find(z_coordinateLn);
                        }
                        else if (line.find(z_coordinate) != string::npos) {
                            x = line.find(z_coordinate);
                        }

                        char last = line[line.length() - 1];


                        string coords;
                        if (line.find("\"") != string::npos && line.find(z_coordinateLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(z_coordinateLn) != string::npos) {

                            size_t y = line.find_first_not_of(" ", x + z_coordinateLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(z_coordinate) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(z_coordinate) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + z_coordinate.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }


                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundAtomsVect[moleculeIteratorRead].size(); j++) {
                                if (compoundAtomsVect[moleculeIteratorRead][j] == currentAtom) {
                                    found = true;
                                    AtomsInfoY[moleculeIteratorRead].push_back(currentAtom);
                                    AtomsInfoY[moleculeIteratorRead].push_back(coords);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundAtomsVect.size(); i++) {
                                for (int j = 0; j < compoundAtomsVect[i].size(); j++) {

                                    if (compoundAtomsVect[i][j] == currentAtom) {
                                        found = true;
                                        AtomsInfoZ[i].push_back(currentAtom);
                                        AtomsInfoZ[i].push_back(coords);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;



                        if (last == '.') {
                            currentAtom = "";
                        }
                    }
                    if ((currentAtom.length() > 2 && line.find(atomName) != string::npos) || (currentAtom.length() > 2 && line.find(atomNameLn) != string::npos)) {

                        size_t x = 0;
                        if (line.find(atomNameLn) != string::npos) {
                            x = line.find(atomNameLn);
                        }
                        else if (line.find(atomName) != string::npos) {
                            x = line.find(atomName);
                        }

                        char last = line[line.length() - 1];


                        string coords;
                        if (line.find("\"") != string::npos && line.find(atomNameLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(atomNameLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + atomNameLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(atomName) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(atomName) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + atomName.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundAtomsVect[moleculeIteratorRead].size(); j++) {
                                if (compoundAtomsVect[moleculeIteratorRead][j] == currentAtom) {
                                    found = true;
                                    AtomsInfoAtom[moleculeIteratorRead].push_back(currentAtom);
                                    AtomsInfoAtom[moleculeIteratorRead].push_back(coords);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundAtomsVect.size(); i++) {
                                for (int j = 0; j < compoundAtomsVect[i].size(); j++) {

                                    if (compoundAtomsVect[i][j] == currentAtom) {
                                        found = true;
                                        AtomsInfoAtom[i].push_back(currentAtom);
                                        AtomsInfoAtom[i].push_back(coords);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;

                        if (last == '.') {
                            currentAtom = "";
                        }
                    }
                    if ((
                        currentAtom.length() == 0
                        && line.find(x_coordinate) != string::npos
                        && line[0] == '_'
                        && line.find(hasAtom) == string::npos
                        && line.find(isCompound) == string::npos
                        && line.find(hasBond) == string::npos
                        && line.find(isBond) == string::npos
                        && line.find(isFirst) == string::npos
                        && line.find(isSecond) == string::npos
                        && line.find(bondsType) == string::npos
                        && line.find(z_coordinate) == string::npos
                        && line.find(y_coordinate) == string::npos
                        && line.find(atomName) == string::npos)
                        ||
                        (currentAtom.length() == 0
                            && line.find(x_coordinateLn) != string::npos
                            && line[0] == '_'
                            && line.find(hasAtomLn) == string::npos
                            && line.find(isCompoundLn) == string::npos
                            && line.find(hasBondLn) == string::npos
                            && line.find(isBondLn) == string::npos
                            && line.find(isFirstLn) == string::npos
                            && line.find(isSecondLn) == string::npos
                            && line.find(bondsTypeLn) == string::npos
                            && line.find(z_coordinateLn) == string::npos
                            && line.find(y_coordinateLn) == string::npos
                            && line.find(atomNameLn) == string::npos)) {



                        size_t y = line.find(" ", 0);
                        currentAtom = line.substr(0, y);

                        size_t x = 0;
                        if (line.find(x_coordinateLn) != string::npos) {
                            x = line.find(x_coordinateLn);
                        }
                        else if (line.find(x_coordinate) != string::npos) {
                            x = line.find(x_coordinate);
                        }

                        char last = line[line.length() - 1];


                        string coords;
                        if (line.find("\"") != string::npos && line.find(x_coordinateLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(x_coordinateLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + x_coordinateLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(x_coordinate) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(x_coordinate) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + x_coordinate.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {

                            for (int j = 0; j < compoundAtomsVect[moleculeIteratorRead].size(); j++) {
                                if (compoundAtomsVect[moleculeIteratorRead][j] == currentAtom) {
                                    found = true;
                                    AtomsInfoX[moleculeIteratorRead].push_back(currentAtom);
                                    AtomsInfoX[moleculeIteratorRead].push_back(coords);
                                    break;
                                }
                            }
                        }
                        else {

                            for (int i = 0; i < compoundAtomsVect.size(); i++) {
                                for (int j = 0; j < compoundAtomsVect[i].size(); j++) {

                                    if (compoundAtomsVect[i][j] == currentAtom) {
                                        found = true;
                                        AtomsInfoX[i].push_back(currentAtom);
                                        AtomsInfoX[i].push_back(coords);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;


                        if (last == '.') {
                            currentAtom = "";
                        }

                    }
                    if ((
                        currentAtom.length() == 0
                        && line.find(y_coordinate) != string::npos
                        && line[0] == '_'
                        && line.find(hasAtom) == string::npos
                        && line.find(isCompound) == string::npos
                        && line.find(hasBond) == string::npos
                        && line.find(isBond) == string::npos
                        && line.find(isFirst) == string::npos
                        && line.find(isSecond) == string::npos
                        && line.find(bondsType) == string::npos
                        && line.find(z_coordinate) == string::npos
                        && line.find(x_coordinate) == string::npos
                        && line.find(atomName) == string::npos)
                        ||
                        (currentAtom.length() == 0
                            && line.find(y_coordinateLn) != string::npos
                            && line[0] == '_'
                            && line.find(hasAtomLn) == string::npos
                            && line.find(isCompoundLn) == string::npos
                            && line.find(hasBondLn) == string::npos
                            && line.find(isBondLn) == string::npos
                            && line.find(isFirstLn) == string::npos
                            && line.find(isSecondLn) == string::npos
                            && line.find(bondsTypeLn) == string::npos
                            && line.find(z_coordinateLn) == string::npos
                            && line.find(x_coordinateLn) == string::npos
                            && line.find(atomNameLn) == string::npos)) {

                        size_t y = line.find(" ", 0);
                        currentAtom = line.substr(0, y);

                        size_t x = 0;
                        if (line.find(y_coordinateLn) != string::npos) {
                            x = line.find(y_coordinateLn);
                        }
                        else if (line.find(y_coordinate) != string::npos) {
                            x = line.find(y_coordinate);
                        }

                        char last = line[line.length() - 1];

                        string coords;
                        if (line.find("\"") != string::npos && line.find(y_coordinateLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(y_coordinateLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + y_coordinateLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(y_coordinate) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(y_coordinate) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + y_coordinate.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundAtomsVect[moleculeIteratorRead].size(); j++) {
                                if (compoundAtomsVect[moleculeIteratorRead][j] == currentAtom) {
                                    found = true;
                                    AtomsInfoY[moleculeIteratorRead].push_back(currentAtom);
                                    AtomsInfoY[moleculeIteratorRead].push_back(coords);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundAtomsVect.size(); i++) {
                                for (int j = 0; j < compoundAtomsVect[i].size(); j++) {

                                    if (compoundAtomsVect[i][j] == currentAtom) {
                                        found = true;
                                        AtomsInfoY[i].push_back(currentAtom);
                                        AtomsInfoY[i].push_back(coords);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;




                        if (last == '.') {
                            currentAtom = "";
                        }

                    }
                    if ((currentAtom.length() == 0 && line.find(z_coordinate) != string::npos && line[0] == '_' && line.find(hasAtom) == string::npos && line.find(isCompound) == string::npos && line.find(hasBond) == string::npos && line.find(isBond) == string::npos && line.find(isFirst) == string::npos && line.find(isSecond) == string::npos && line.find(bondsType) == string::npos && line.find(x_coordinate) == string::npos && line.find(y_coordinate) == string::npos && line.find(atomName) == string::npos) || (currentAtom.length() == 0 && line.find(z_coordinateLn) != string::npos && line[0] == '_' && line.find(hasAtomLn) == string::npos && line.find(isCompoundLn) == string::npos && line.find(hasBondLn) == string::npos && line.find(isBondLn) == string::npos && line.find(isFirstLn) == string::npos && line.find(isSecondLn) == string::npos && line.find(bondsTypeLn) == string::npos && line.find(x_coordinateLn) == string::npos && line.find(y_coordinateLn) == string::npos && line.find(atomNameLn) == string::npos)) {

                        size_t y = line.find(" ", 0);
                        currentAtom = line.substr(0, y);

                        size_t x = 0;
                        if (line.find(z_coordinateLn) != string::npos) {
                            x = line.find(z_coordinateLn);
                        }
                        else if (line.find(z_coordinate) != string::npos) {
                            x = line.find(z_coordinate);
                        }

                        char last = line[line.length() - 1];


                        string coords;
                        if (line.find("\"") != string::npos && line.find(z_coordinateLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(z_coordinateLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + z_coordinateLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(z_coordinate) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(z_coordinate) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + z_coordinate.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundAtomsVect[moleculeIteratorRead].size(); j++) {
                                if (compoundAtomsVect[moleculeIteratorRead][j] == currentAtom) {
                                    found = true;
                                    AtomsInfoZ[moleculeIteratorRead].push_back(currentAtom);
                                    AtomsInfoZ[moleculeIteratorRead].push_back(coords);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundAtomsVect.size(); i++) {
                                for (int j = 0; j < compoundAtomsVect[i].size(); j++) {

                                    if (compoundAtomsVect[i][j] == currentAtom) {
                                        found = true;
                                        AtomsInfoZ[i].push_back(currentAtom);
                                        AtomsInfoZ[i].push_back(coords);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }
                        found = false;




                        if (last == '.') {
                            currentAtom = "";
                        }

                    }
                    if ((currentAtom.length() == 0 && line.find(atomName) != string::npos && line[0] == '_' && line.find(hasAtom) == string::npos && line.find(isCompound) == string::npos && line.find(hasBond) == string::npos && line.find(isBond) == string::npos && line.find(isFirst) == string::npos && line.find(isSecond) == string::npos && line.find(bondsType) == string::npos && line.find(x_coordinate) == string::npos && line.find(y_coordinate) == string::npos && line.find(z_coordinate) == string::npos) || (currentAtom.length() == 0 && line.find(atomNameLn) != string::npos && line[0] == '_' && line.find(hasAtomLn) == string::npos && line.find(isCompoundLn) == string::npos && line.find(hasBondLn) == string::npos && line.find(isBondLn) == string::npos && line.find(isFirstLn) == string::npos && line.find(isSecondLn) == string::npos && line.find(bondsTypeLn) == string::npos && line.find(x_coordinateLn) == string::npos && line.find(z_coordinateLn) == string::npos && line.find(y_coordinateLn) == string::npos)) {

                        size_t y = line.find(" ", 0);
                        currentAtom = line.substr(0, y);

                        size_t x = 0;
                        if (line.find(atomNameLn) != string::npos) {
                            x = line.find(atomNameLn);
                        }
                        else if (line.find(atomName) != string::npos) {
                            x = line.find(atomName);
                        }

                        char last = line[line.length() - 1];


                        string coords;
                        if (line.find("\"") != string::npos && line.find(atomNameLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(atomNameLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + atomNameLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(atomName) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            coords = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(atomName) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + atomName.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            coords = line.substr(y, (u - y));


                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundAtomsVect[moleculeIteratorRead].size(); j++) {
                                if (compoundAtomsVect[moleculeIteratorRead][j] == currentAtom) {
                                    found = true;
                                    AtomsInfoAtom[moleculeIteratorRead].push_back(currentAtom);
                                    AtomsInfoAtom[moleculeIteratorRead].push_back(coords);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundAtomsVect.size(); i++) {
                                for (int j = 0; j < compoundAtomsVect[i].size(); j++) {

                                    if (compoundAtomsVect[i][j] == currentAtom) {
                                        found = true;
                                        AtomsInfoAtom[i].push_back(currentAtom);
                                        AtomsInfoAtom[i].push_back(coords);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;


                        if (last == '.') {
                            currentAtom = "";
                        }

                    }
                    ////////------------------Retrieving Bonds INFO -----------------------//////
                     //Bonds INFOS
                    if ((line[0] == '_' && line.find(isBond) != string::npos && line[line.length() - 1] == ';') || (line[0] == '_' && line.find(isBondLn) != string::npos && line[line.length() - 1] == ';') || (line[0] == '_' && line.find(isBond) != string::npos && line[line.length() - 1] == '.') || (line[0] == '_' && line.find(isBondLn) != string::npos && line[line.length() - 1] == '.')) {

                        size_t x = line.find(" ", 0);

                        currentBond = line.substr(0, x);
                        char last = line[line.length() - 1];
                        if (last == '.') {
                            currentBond = "";
                        }
                    }

                    if ((currentBond.length() > 2 && line.find(isFirst) != string::npos) || (currentBond.length() > 2 && line.find(isFirstLn) != string::npos)) {

                        size_t x = 0;
                        if (line.find(isFirstLn) != string::npos) {
                            x = line.find(isFirstLn);
                        }
                        else if (line.find(isFirst) != string::npos) {
                            x = line.find(isFirst);
                        }

                        char last = line[line.length() - 1];



                        string first;
                        if (line.find("\"") != string::npos && line.find(isFirstLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            first = line.substr(y + 1, (u - y) - 1);


                        }
                        else if (line.find(isFirstLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + isFirstLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            first = line.substr(y, (u - y));



                        }
                        if (line.find("\"") != string::npos && line.find(isFirst) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            first = line.substr(y + 1, (u - y) - 1);


                        }
                        else if (line.find(isFirst) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + isFirst.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            first = line.substr(y, (u - y));



                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundBondsVect[moleculeIteratorRead].size(); j++) {

                                if (compoundBondsVect[moleculeIteratorRead][j] == currentBond) {
                                    found = true;
                                    BondsInfoFirst[moleculeIteratorRead].push_back(currentBond);
                                    BondsInfoFirst[moleculeIteratorRead].push_back(first);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundBondsVect.size(); i++) {
                                for (int j = 0; j < compoundBondsVect[i].size(); j++) {

                                    if (compoundBondsVect[i][j] == currentBond) {
                                        found = true;
                                        BondsInfoFirst[i].push_back(currentBond);
                                        BondsInfoFirst[i].push_back(first);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }
                        found = false;


                        if (last == '.') {
                            currentBond = "";
                        }
                    }
                    if ((currentBond.length() > 2 && line.find(isSecond) != string::npos) || (currentBond.length() > 2 && line.find(isSecondLn) != string::npos)) {

                        size_t x = 0;
                        if (line.find(isSecondLn) != string::npos) {
                            x = line.find(isSecondLn);
                        }
                        else if (line.find(isSecond) != string::npos) {
                            x = line.find(isSecond);
                        }

                        char last = line[line.length() - 1];


                        string second;
                        if (line.find("\"") != string::npos && line.find(isSecondLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            second = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(isSecondLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + isSecondLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            second = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(isSecond) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            second = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(isSecond) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + isSecond.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            second = line.substr(y, (u - y));


                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundBondsVect[moleculeIteratorRead].size(); j++) {

                                if (compoundBondsVect[moleculeIteratorRead][j] == currentBond) {
                                    found = true;
                                    BondsInfoSecond[moleculeIteratorRead].push_back(currentBond);
                                    BondsInfoSecond[moleculeIteratorRead].push_back(second);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundBondsVect.size(); i++) {
                                for (int j = 0; j < compoundBondsVect[i].size(); j++) {

                                    if (compoundBondsVect[i][j] == currentBond) {
                                        found = true;
                                        BondsInfoSecond[i].push_back(currentBond);
                                        BondsInfoSecond[i].push_back(second);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;


                        if (last == '.') {
                            currentBond = "";
                        }
                    }
                    if ((currentBond.length() > 2 && line.find(bondsType) != string::npos) || (currentBond.length() > 2 && line.find(bondsTypeLn) != string::npos)) {

                        size_t x = 0;
                        if (line.find(bondsTypeLn) != string::npos) {
                            x = line.find(bondsTypeLn);
                        }
                        else if (line.find(bondsType) != string::npos) {
                            x = line.find(bondsType);
                        }

                        char last = line[line.length() - 1];


                        string type;
                        if (line.find("\"") != string::npos && line.find(bondsTypeLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            type = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(bondsTypeLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + bondsTypeLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            type = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(bondsType) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            type = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(bondsType) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + bondsType.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            type = line.substr(y, (u - y));


                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundBondsVect[moleculeIteratorRead].size(); j++) {

                                if (compoundBondsVect[moleculeIteratorRead][j] == currentBond) {
                                    found = true;
                                    BondsInfoType[moleculeIteratorRead].push_back(currentBond);
                                    BondsInfoType[moleculeIteratorRead].push_back(type);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundBondsVect.size(); i++) {
                                for (int j = 0; j < compoundBondsVect[i].size(); j++) {

                                    if (compoundBondsVect[i][j] == currentBond) {
                                        found = true;
                                        BondsInfoType[i].push_back(currentBond);
                                        BondsInfoType[i].push_back(type);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;


                        if (last == '.') {
                            currentBond = "";
                        }
                    }
                    if ((currentBond.length() == 0
                        && line.find(isFirst) != string::npos && line[0] == '_' && line.find(hasAtom) == string::npos && line.find(isCompound) == string::npos && line.find(hasBond) == string::npos && line.find(isBond) == string::npos && line.find(isAtom) == string::npos && line.find(x_coordinate) == string::npos && line.find(y_coordinate) == string::npos && line.find(z_coordinate) == string::npos && line.find(atomName) == string::npos && line.find(isSecond) == string::npos && line.find(bondsType) == string::npos) || (currentBond.length() == 0 && line.find(isFirstLn) != string::npos && line[0] == '_' && line.find(hasAtomLn) == string::npos && line.find(isCompoundLn) == string::npos && line.find(hasBondLn) == string::npos && line.find(isBondLn) == string::npos && line.find(isAtomLn) == string::npos && line.find(x_coordinateLn) == string::npos && line.find(y_coordinateLn) == string::npos && line.find(z_coordinateLn) == string::npos && line.find(atomNameLn) == string::npos && line.find(isSecondLn) == string::npos && line.find(bondsTypeLn) == string::npos)) {

                        size_t y = line.find(" ", 0);
                        currentBond = line.substr(0, y);

                        size_t x = 0;
                        if (line.find(isFirstLn) != string::npos) {
                            x = line.find(isFirstLn);
                        }
                        else if (line.find(isFirst) != string::npos) {
                            x = line.find(isFirst);
                        }

                        char last = line[line.length() - 1];


                        string first;
                        if (line.find("\"") != string::npos && line.find(isFirstLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            first = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(isFirstLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + isFirstLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            first = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(isFirst) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            first = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(isFirst) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + isFirst.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            first = line.substr(y, (u - y));

                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundBondsVect[moleculeIteratorRead].size(); j++) {

                                if (compoundBondsVect[moleculeIteratorRead][j] == currentBond) {
                                    found = true;
                                    BondsInfoFirst[moleculeIteratorRead].push_back(currentBond);
                                    BondsInfoFirst[moleculeIteratorRead].push_back(first);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundBondsVect.size(); i++) {
                                for (int j = 0; j < compoundBondsVect[i].size(); j++) {

                                    if (compoundBondsVect[i][j] == currentBond) {
                                        found = true;
                                        BondsInfoFirst[i].push_back(currentBond);
                                        BondsInfoFirst[i].push_back(first);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;


                        if (last == '.') {
                            currentBond = "";
                        }
                    }
                    if ((currentBond.length() == 0 && line.find(isSecond) != string::npos && line[0] == '_' && line.find(hasAtom) == string::npos && line.find(isCompound) == string::npos && line.find(hasBond) == string::npos && line.find(isBond) == string::npos && line.find(isAtom) == string::npos && line.find(x_coordinate) == string::npos && line.find(y_coordinate) == string::npos && line.find(z_coordinate) == string::npos && line.find(atomName) == string::npos && line.find(isFirst) == string::npos && line.find(bondsType) == string::npos) || (currentBond.length() == 0 && line.find(isSecondLn) != string::npos && line[0] == '_' && line.find(hasAtomLn) == string::npos && line.find(isCompoundLn) == string::npos && line.find(hasBondLn) == string::npos && line.find(isBondLn) == string::npos && line.find(isAtomLn) == string::npos && line.find(x_coordinateLn) == string::npos && line.find(y_coordinateLn) == string::npos && line.find(z_coordinateLn) == string::npos && line.find(atomNameLn) == string::npos && line.find(isFirstLn) == string::npos && line.find(bondsTypeLn) == string::npos)) {

                        size_t y = line.find(" ", 0);
                        currentBond = line.substr(0, y);

                        size_t x = 0;
                        if (line.find(isSecondLn) != string::npos) {
                            x = line.find(isSecondLn);
                        }
                        else if (line.find(isSecond) != string::npos) {
                            x = line.find(isSecond);
                        }

                        char last = line[line.length() - 1];

                        string second;
                        if (line.find("\"") != string::npos && line.find(isSecondLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            second = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(isSecondLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + isSecondLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            second = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(isSecond) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            second = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(isSecond) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + isSecond.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            second = line.substr(y, (u - y));

                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundBondsVect[moleculeIteratorRead].size(); j++) {

                                if (compoundBondsVect[moleculeIteratorRead][j] == currentBond) {
                                    found = true;
                                    BondsInfoSecond[moleculeIteratorRead].push_back(currentBond);
                                    BondsInfoSecond[moleculeIteratorRead].push_back(second);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundBondsVect.size(); i++) {
                                for (int j = 0; j < compoundBondsVect[i].size(); j++) {

                                    if (compoundBondsVect[i][j] == currentBond) {
                                        found = true;
                                        BondsInfoSecond[i].push_back(currentBond);
                                        BondsInfoSecond[i].push_back(second);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;


                        if (last == '.') {
                            currentBond = "";
                        }

                    }
                    if ((currentBond.length() == 0 && line.find(bondsType) != string::npos && line[0] == '_' && line.find(hasAtom) == string::npos && line.find(isCompound) == string::npos && line.find(hasBond) == string::npos && line.find(isBond) == string::npos && line.find(isAtom) == string::npos && line.find(x_coordinate) == string::npos && line.find(y_coordinate) == string::npos && line.find(z_coordinate) == string::npos && line.find(atomName) == string::npos && line.find(isFirst) == string::npos && line.find(isSecond) == string::npos) || (currentBond.length() == 0 && line.find(bondsTypeLn) != string::npos && line[0] == '_' && line.find(hasAtomLn) == string::npos && line.find(isCompoundLn) == string::npos && line.find(hasBondLn) == string::npos && line.find(isBondLn) == string::npos && line.find(isAtomLn) == string::npos && line.find(x_coordinateLn) == string::npos && line.find(y_coordinateLn) == string::npos && line.find(z_coordinateLn) == string::npos && line.find(atomNameLn) == string::npos && line.find(isFirstLn) == string::npos && line.find(isSecondLn) == string::npos)) {

                        size_t y = line.find(" ", 0);
                        currentBond = line.substr(0, y);

                        size_t x = 0;
                        if (line.find(bondsTypeLn) != string::npos) {
                            x = line.find(bondsTypeLn);
                        }
                        else if (line.find(bondsType) != string::npos) {
                            x = line.find(bondsType);
                        }

                        char last = line[line.length() - 1];

                        string type;
                        if (line.find("\"") != string::npos && line.find(bondsTypeLn) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            type = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(bondsTypeLn) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + bondsTypeLn.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            type = line.substr(y, (u - y));


                        }
                        if (line.find("\"") != string::npos && line.find(bondsType) != string::npos) {

                            size_t y = line.find("\"", x + 1);
                            size_t u = line.find("\"", y + 1);
                            type = line.substr(y + 1, (u - y) - 1);

                        }
                        else if (line.find(bondsType) != string::npos) {

                            size_t y = line.find_first_not_of(" \t", x + bondsType.size() + 1);
                            size_t u = line.find(" ", y + 1);

                            type = line.substr(y, (u - y));


                        }

                        if (pConv->IsOption("S", OBConversion::INOPTIONS))
                        {
                            for (int j = 0; j < compoundBondsVect[moleculeIteratorRead].size(); j++) {

                                if (compoundBondsVect[moleculeIteratorRead][j] == currentBond) {
                                    found = true;
                                    BondsInfoType[moleculeIteratorRead].push_back(currentBond);
                                    BondsInfoType[moleculeIteratorRead].push_back(type);
                                    break;
                                }
                            }
                        }
                        else {
                            for (int i = 0; i < compoundBondsVect.size(); i++) {
                                for (int j = 0; j < compoundBondsVect[i].size(); j++) {

                                    if (compoundBondsVect[i][j] == currentBond) {
                                        found = true;
                                        BondsInfoType[i].push_back(currentBond);
                                        BondsInfoType[i].push_back(type);
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }

                        found = false;


                        if (last == '.') {
                            currentBond = "";
                        }


                    }

                }

                /* For multi-molecule formats, leave the input stream at the start of the
                   next molecule, ready for this routine to be called again.

                /* Return true if ok. Returning false means discard the OBMol and stop
                   converting, unless the -e option is set. With a multi-molecule inputstream
                   this will skip the current molecule and continue with the next, if SkipObjects()
                   has been defined. If it has not, and continuation after errors is still required,
                   it is necessary to leave the input stream at the beginning of next object when
                   returning false;*/

            }
            cerr << "Finished processing file" << "\n";
            isFirstRead = false;
        }

        atmNmIndx = 1;
        int goo = 1;
        for (int j = 0; j < compoundAtomsIndxVect[moleculeIterator].size(); j++) {
            if (j == goo) {
                compoundAtomsIndxVect[moleculeIterator][j] = std::to_string(atmNmIndx);
                atmNmIndx++;
                goo += 2;
            }
        }

        unsigned int begin, end, order, flag;

        pmol->BeginModify();
        //
        // ATOM Block
        //

        int atomsReserve = compoundAtomsVect[moleculeIterator].size();
        mol.ReserveAtoms(atomsReserve);

        double x, y, z;
        string symbol;
        int massdiff, charge, stereo;

        for (int f = 0; f < compoundAtomsVect[moleculeIterator].size(); f++) {

            massdiff = charge = 0;

            string currAtm = compoundAtomsVect[moleculeIterator][f];
            for (int g = 1; g < AtomsInfoX[moleculeIterator].size(); g++) {

                if (AtomsInfoX[moleculeIterator][g - 1] == currAtm) {

                    x = stod(AtomsInfoX[moleculeIterator][g]);
                    foundX = true;
                    break;
                }
            }
            for (int h = 1; h < AtomsInfoY[moleculeIterator].size(); h++) {
                if (AtomsInfoY[moleculeIterator][h - 1] == currAtm) {

                    y = stod(AtomsInfoY[moleculeIterator][h]);
                    foundY = true;
                    break;
                }
            }
            for (int j = 1; j < AtomsInfoZ[moleculeIterator].size(); j++) {
                if (AtomsInfoZ[moleculeIterator][j - 1] == currAtm) {

                    z = stod(AtomsInfoZ[moleculeIterator][j]);
                    foundZ = true;
                    break;
                }
            }
            for (int k = 1; k < AtomsInfoAtom[moleculeIterator].size(); k++) {
                if (AtomsInfoAtom[moleculeIterator][k - 1] == currAtm) {

                    symbol = AtomsInfoAtom[moleculeIterator][k];
                    foundAtm = true;
                    break;
                }
            }

            if (foundX && foundY && foundZ && foundAtm) {
                OBAtom* patom = mol.NewAtom();
                patom->SetVector(x, y, z);
                patom->SetAtomicNum(OBElements::GetAtomicNum(symbol.c_str()));
                foundX = false;
                foundY = false;
                foundZ = false;
                foundAtm = false;
            }

        }

        //
        // Bond Block
        //
        foundBndFir = false;
        foundBndSec = false;
        foundBndType = false;

        stereo = 0;
        flag = 0;
        int iterG = 1;

        for (int g = 1; g < BondsInfoFirst[moleculeIterator].size(); g++) {

            if (g == iterG) {
                int iterH = 1;
                for (int h = 1; h < compoundAtomsIndxVect[moleculeIterator].size(); h++) {

                    if (h == iterH) {
                        iterH += 2;
                        if (BondsInfoFirst[moleculeIterator][g] == compoundAtomsIndxVect[moleculeIterator][h - 1]) {

                            string convBegin = compoundAtomsIndxVect[moleculeIterator][h];
                            begin = stoi(convBegin);
                            foundBndFir = true;
                        }
                        if (BondsInfoSecond[moleculeIterator][g] == compoundAtomsIndxVect[moleculeIterator][h - 1]) {

                            string convEnd = compoundAtomsIndxVect[moleculeIterator][h];
                            end = stoi(convEnd);
                            foundBndSec = true;
                        }

                        if (BondsInfoType[moleculeIterator].size() >= h) {
                            if (BondsInfoFirst[moleculeIterator][g - 1] == BondsInfoType[moleculeIterator][h - 1]) {

                                if (BondsInfoType[moleculeIterator][h] == ":single") {

                                    order = 1;
                                }
                                else if (BondsInfoType[moleculeIterator][h] == ":double") {
                                    order = 2;

                                }
                                else if (BondsInfoType[moleculeIterator][h] == ":aromatic") {
                                    order = 4;

                                }
                                else if (BondsInfoType[moleculeIterator][h] == "<https://ii.uwb.edu.pl/ctld#single>") {
                                    order = 1;

                                }
                                else if (BondsInfoType[moleculeIterator][h] == "<https://ii.uwb.edu.pl/ctld#double>") {
                                    order = 2;

                                }
                                else if (BondsInfoType[moleculeIterator][h] == "<https://ii.uwb.edu.pl/ctld#aromatic>") {
                                    order = 4;

                                }
                                else {
                                    order = 1;

                                }
                                foundBndType = true;
                            }
                        }
                        if (foundBndType && foundBndSec && foundBndFir) {
                            if (!mol.AddBond(begin, end, order)) {
                                errorMsg << "WARNING: Problems reading a CTLD file\n";
                                errorMsg << "Invalid bond specification\n";
                                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
                                return false;
                            }
                            foundBndType = false;
                            foundBndSec = false;
                            foundBndFir = false;
                            break;
                        }
                    }
                }
                iterG += 2;
            }

        }

        /** Parse the input stream and use the OpenBabel API to populate the OBMol **/

        // To use an input option
        // if(pConv->IsOption("s",OBConversion::INOPTIONS))
        // {
        // 	//Code for when -as is specified
        // }

        /* If the molecule has other than 3D coordinates for its atoms, it
        is necessary to set the dimension to 0, or 2 */
        // int dim;
        pmol->SetDimension(2);

        pmol->EndModify();

        moleculeIterator += 1;

        if (moleculeIterator < compoundWithIndex.size()) {

            ifs.clear();
            ifs.seekg(0);

        }
        if (moleculeIterator == compoundWithIndex.size()) {
            ifs.seekg(EOF);
        }

        return true;
    }

    ////////////////////////////////////////////////////////////////

    bool CTLDFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
    {
        OBMol* pmol = dynamic_cast<OBMol*>(pOb);

        if (pmol == nullptr)
            return false;
        //Define some references so we can use the old parameter names
        ostream& ofs = *pConv->GetOutStream();
        OBMol& mol = *pmol;

        if (mol.GetDimension() == 0)
        {
            obErrorLog.ThrowError(__FUNCTION__, "No 2D or 3D coordinates exist. Stereochemical information will", obWarning, onceOnly);
        }

        string dimension("2D");
        if (mol.GetDimension() == 3)
            dimension = "3D";

        char buff[BUFF_SIZE];

        if (mol.NumAtoms() > 999 || mol.NumBonds() > 999) { // Three digits!
            stringstream errorMsg;
            errorMsg << "CTLD file conversion failed: Molecule is too large to convert." << endl;
            errorMsg << "  File format is limited to 999 atoms or bonds." << endl;
            errorMsg << "  Molecule size: " << mol.NumAtoms() << " atoms ";
            errorMsg << "and " << mol.NumBonds() << " bonds." << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
            return false;
        }

        set<OBBond*> unspec_ctstereo = GetUnspecifiedCisTrans(mol);
        map<OBAtom*, Parity> parity;
        map<OBBond*, OBStereo::Ref> from;
        map<OBBond*, OBStereo::Ref>::const_iterator from_cit;
        GetParity(mol, parity);
        if (!prefixApplied) {
            if (writePrefix != "") {
                ofs << writePrefix << endl;
                prefixApplied = true;
            }
        }
        writeMolIterator++;
        int currentMol = writeMolIterator;
        OBConversion conv;
        if (pConv->IsOption("E", OBConversion::OUTOPTIONS))
        {
            OBFormat* canSMIFormat = conv.FindFormat("can");
            OBFormat* inchiFormat = conv.FindFormat("inchi");
            OBFormat* inchiKeyFormat = conv.FindFormat("inchikey");
            snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
            ofs << buff << " " << rdfType << " " << writeCompound << writeEndPoint << endl;
            ofs << buff << " " << enrichFormula << " \"" << mol.GetFormula() << "\"" << writeXMLstring << writeEndPoint << endl;
            ofs << buff << " " << enrichMolecularWeight << " " << writeColon << mol.GetMolWt() << writeColon << writeXMLdecimal << writeEndPoint << endl;
            ofs << buff << " " << enrichExactMass << " " << writeColon << mol.GetExactMass() << writeColon << writeXMLdecimal << writeEndPoint << endl;
            OBDescriptor* pDesc;
            pDesc = OBDescriptor::FindType("logP");
            if (pDesc)
                ofs << buff << " " << enrichLogP << " " << writeColon << pDesc->Predict(&mol) << writeColon << writeXMLdecimal << writeEndPoint << endl;
            pDesc = OBDescriptor::FindType("TPSA");
            if (pDesc)
                ofs << buff << " " << enrichTpsa << " " << writeColon << pDesc->Predict(&mol) << writeColon << writeXMLdecimal << writeEndPoint << endl;
            pDesc = OBDescriptor::FindType("MR");
            if (pDesc)
                ofs << buff << " " << enrichMolecularRefractivity << " " << writeColon << pDesc->Predict(&mol) << writeColon << writeXMLdecimal << writeEndPoint << endl;
            string smilesString = "-";

            if (canSMIFormat) {
                conv.SetOutFormat(canSMIFormat);

                smilesString = conv.WriteString(&mol);

                if (smilesString.length() == 0) {
                    smilesString = "-";
                }
                else {
                     EscapeRDF(smilesString);
                    ofs << buff << " " << enrichSmiles << " \"" << smilesString << "\"" << writeXMLstring << writeEndPoint << endl;
                }
            }
            string inchiString = "-";
            if (inchiFormat) {
                conv.SetOutFormat(inchiFormat);
                inchiString = conv.WriteString(&mol);
                if (inchiString.length() == 0) {
                    inchiString = "-";
                }
                else {
                    EscapeRDF(inchiString);
                    ofs << buff << " " << enrichInchi << " \"" << inchiString << "\"" << writeXMLstring << writeEndPoint << endl;
                }
            }
            string inchiKeyString = "-";
            if (inchiFormat) {
                conv.SetOutFormat(inchiKeyFormat);
                inchiKeyString = conv.WriteString(&mol);
                if (inchiKeyString.length() == 0) {
                    inchiKeyString = "-";
                }
                else {
                    EscapeRDF(inchiKeyString);
                    ofs << buff << " "  << enrichInchiKey << " \"" << inchiKeyString << "\"" << writeXMLstring << writeEndPoint << endl;
                }
            }
            ofs << writeline;
        }
        else {
            snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
            ofs << buff << " " << rdfType << " " << writeCompound << writeEndPoint << endl;
            ofs << writeline;
        }
        OBAtom *atom;
        vector<OBAtom*>::iterator i;
        unsigned int aclass = 0;
        int charge = 0;
        for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i)) {

            Parity stereo = NotStereo;
            if (parity.find(atom) != parity.end())
                stereo = parity[atom];

            writeMolIterator++;

            snprintf(buff, BUFF_SIZE, "_:b%i", currentMol);
            ofs << buff << " " << writeHasAtom << " ";
            snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
            ofs << buff << writeEndPoint << endl;
            ofs << writeline;

            snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
            ofs << buff << " " << rdfType << " " << writeAtom << writeEndPoint << endl;

            snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
            ofs << buff << " " << writeatom << " ";
            snprintf(buff, BUFF_SIZE, "\"%s\"", OBElements::GetSymbol(atom->GetAtomicNum()));
            ofs << buff << writeXMLstring << writeEndPoint << endl;

            snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
            ofs << buff << " " << writexCoordinate << " ";
            snprintf(buff, BUFF_SIZE, "%.4f", atom->GetX());
            ofs << writeColon << buff << writeColon << writeXMLdecimal << writeEndPoint << endl;

            snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
            ofs << buff << " " << writeyCoordinate << " ";
            snprintf(buff, BUFF_SIZE, "%.4f", atom->GetY());
            ofs << writeColon << buff << writeColon << writeXMLdecimal << writeEndPoint << endl;

            snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
            ofs << buff << " " << writezCoordinate << " ";
            snprintf(buff, BUFF_SIZE, "%.4f", atom->GetZ());
            ofs << writeColon << buff << writeColon << writeXMLdecimal << writeEndPoint << endl;

            if (stereo != 0) {
                snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
                ofs << buff << " " << writeStereo << " " << writeColon << stereo << writeColon << writeXMLdecimal << writeEndPoint << endl;
            }
            ofs << writeline;


        }
        OBAtom *nbr;
        OBBond *bond;
        vector<OBBond*>::iterator j;
        int bondline = 0;
        vector<int> zbos;
        for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i)) {
            for (nbr = atom->BeginNbrAtom(j); nbr; nbr = atom->NextNbrAtom(j)) {
                bond = (OBBond*)*j;
                from_cit = from.find(bond);
                // If the bond has *calculated* stereodirectionality, ensure that the start point
                // is at the 'from' atom. Otherwise, just ensure that the start atom
                // is the 'begin atom' of the bond (so that stereodirectionality that was
                // read in [rather than calculated] will be correct).
                if ((from_cit == from.end() && atom->GetIdx() == bond->GetBeginAtomIdx()) || (from_cit != from.end() && from_cit->second == atom->GetId())) {

                    int stereo = 0;
                    if (mol.GetDimension() == 2 && pConv->IsOption("w", pConv->OUTOPTIONS) != nullptr) {
                        if (bond->IsWedge())
                            stereo = 1;
                        else if (bond->IsHash())
                            stereo = 6;
                        else if (bond->IsWedgeOrHash())
                            stereo = 4;
                    }

                    if (unspec_ctstereo.find(bond) != unspec_ctstereo.end())
                        stereo = 3;


                    int firstAt = atom->GetIdx(); // begin atom number
                    int secndAt = nbr->GetIdx(); // end atom number
                    int firstAtMod = currentMol + firstAt;
                    int secndAtMod = currentMol + secndAt;
                    int bondTypeNum = (int)(bond->GetBondOrder());
                    string bondTypeMod = writeSingle;
                    switch (bondTypeNum) {
                    case 2:
                        bondTypeMod = writeDouble;
                        break;
                    default:
                        bondTypeMod = writeSingle;
                    }

                    writeMolIterator++;

                    snprintf(buff, BUFF_SIZE, "_:b%i", currentMol);
                    ofs << buff << " " << writeHasBond << " ";
                    snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
                    ofs << buff << writeEndPoint << endl;
                    ofs << writeline;

                    snprintf(buff, BUFF_SIZE, "_:b%i",
                        writeMolIterator);
                    ofs << buff << " " << rdfType << " " << writeBond << writeEndPoint << endl;

                    snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
                    ofs << buff << " " << writeFirstAtom << " ";
                    snprintf(buff, BUFF_SIZE, "_:b%i", firstAtMod);
                    ofs << buff << writeEndPoint << endl;

                    snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
                    ofs << buff << " " << writeSecondAtom << " ";
                    snprintf(buff, BUFF_SIZE, "_:b%i", secndAtMod);
                    ofs << buff << writeEndPoint << endl;


                    snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
                    ofs << buff << " " << writeType << " " << bondTypeMod << writeEndPoint << endl;

                    if (stereo != 0) {
                        snprintf(buff, BUFF_SIZE, "_:b%i", writeMolIterator);
                        ofs << buff << " " << writeStereo << " " << writeColon << stereo << writeColon << writeEndPoint << endl;
                    }
                    ofs << writeline;

                }
            }
        }

        return true; //or false to stop converting
    }
    void CTLDFormat::EscapeRDF(string& str) {
        size_t start_pos = 0;
        while ((start_pos = str.find("\\", start_pos)) != string::npos) {
            str.replace(start_pos, 1, "\\\\");
            start_pos += 2;
        }
        if (!str.empty() && str[str.length() - 1] == '\n') {
            str.erase(str.length() - 1);
        }
        if (!str.empty() && str[str.length() - 1] == '\t') {
            str.erase(str.length() - 1);
        }
        if (!str.empty() && str[str.length() - 1] == '\r') {
            str.erase(str.length() - 1);
        }
        if (!str.empty() && str[str.length() - 1] == '\t') {
            str.erase(str.length() - 1);
        }
        
    }
    void CTLDFormat::GetParity(OBMol& mol, map<OBAtom*, CTLDFormat::Parity>& parity)
    {
        // This loop sets the atom parity for each tet center
        std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
        for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
            if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
                OBTetrahedralStereo* ts = dynamic_cast<OBTetrahedralStereo*>(*data);

                OBTetrahedralStereo::Config cfg = ts->GetConfig();

                Parity atomparity = Unknown;
                if (cfg.specified && cfg.winding != OBStereo::UnknownWinding) {
                    // If, when looking towards the maxref, the remaining refs increase in number
                    // clockwise, parity is 1 (Parity::Clockwise). Note that Implicit Refs and Hydrogens
                    // should be treated considered the maxref if present.
                    OBStereo::Refs refs = cfg.refs;

                    unsigned long maxref = OBStereo::NoRef;
                    // Search for an explicit Hydrogen in the cfg refs...
                    if (cfg.from != OBStereo::ImplicitRef && mol.GetAtomById(cfg.from)->GetAtomicNum() == OBElements::Hydrogen)
                        maxref = cfg.from;
                    else
                        for (OBStereo::RefIter ref_it = refs.begin(); ref_it != refs.end(); ++ref_it)
                            if ((*ref_it) != OBStereo::ImplicitRef && mol.GetAtomById(*ref_it)->GetAtomicNum() == OBElements::Hydrogen)
                                maxref = *ref_it;
                    // ...otherwise, find the maximum ref (note that ImplicitRef will be max if present)
                    if (maxref == OBStereo::NoRef)
                        maxref = std::max(*(std::max_element(refs.begin(), refs.end())), cfg.from);

                    // Get a new cfg and refs looking towards the maxref
                    cfg = ts->GetConfig(maxref, OBStereo::Clockwise, OBStereo::ViewTowards);
                    int inversions = OBStereo::NumInversions(cfg.refs);

                    // If they were in increasing order, inversions would be 0 or some even value
                    if (inversions % 2 == 0)
                        atomparity = Clockwise;
                    else
                        atomparity = AntiClockwise;
                }
                parity[mol.GetAtomById(cfg.center)] = atomparity;
            }
    }
}
//namespace OpenBabel
