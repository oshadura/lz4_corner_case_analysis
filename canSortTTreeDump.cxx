//
// Usage examples:
// root [0] .x canSortTTreeDump.cxx
// root [0] .x canSortTTreeDump.cxx++
// root [0] .x canSortTTreeDump.cxx("test.root")
// root [0] .x canSortTTreeDump.cxx++("test.root")
// root [0] .x canSortTTreeDump.cxx("test.root", "canSort")
// root [0] .x canSortTTreeDump.cxx++("test.root", "canSort")
//

#include "TFile.h"
#include "TTree.h"

#include <iostream>

void canSortTTreeDump(const char *filename = "test.root",
                      const char *treename = "canSort")
{
  if (!(filename && filename[0])) return; // just a precaution
  if (!(treename && treename[0])) return; // just a precaution
  
  TFile *f = TFile::Open(filename);
  if ((!f) || f->IsZombie()) { delete f; return; } // just a precaution
  
  TTree *t; f->GetObject(treename, t);
  if (!t) { delete f; return; } // just a precaution
  t->Print();
  
  const Int_t _MAX_No_Clover_ = 16;
  Int_t No_Clover;
#if 1 /* 0 or 1 */
  Int_t clov_mult, Clover[(_MAX_No_Clover_)];
  ULong64_t CLT[(_MAX_No_Clover_)], Eclab[(_MAX_No_Clover_)];
#elif 1 /* 0 or 1 */
  Short_t clov_mult, Clover[(_MAX_No_Clover_)];
  UShort_t CLT[(_MAX_No_Clover_)], Eclab[(_MAX_No_Clover_)];
#else /* 0 or 1 */
  Char_t clov_mult, Clover[(_MAX_No_Clover_)];
  UShort_t CLT[(_MAX_No_Clover_)], Eclab[(_MAX_No_Clover_)];
#endif /* 0 or 1 */
  
  t->SetBranchAddress("clov_mult", &clov_mult);
  t->SetBranchAddress("No_Clover", &No_Clover);
  t->SetBranchAddress("Clover", Clover);
  t->SetBranchAddress("CLT", CLT);
  t->SetBranchAddress("Eclab", Eclab);
  
  Long64_t n = t->GetEntries();
#if 1 /* 0 or 1 */
  if (n > 33) n = 33; // process the first 33 entries only
#endif /* 0 or 1 */
  for (Long64_t i = 0; i < n; i++) {
    if (t->GetEntry(i) < 1) break; // "break" if any problem met
    if (No_Clover > (_MAX_No_Clover_)) { // just a precaution
      std::cout << "We are doomed! ... No_Clover = " << No_Clover << std::endl;
      delete f; // automatically deletes "t", too
      return;
    }
#if 1 /* 0 or 1 */
    // dump all entries
    std::cout << i << " : "
              << ((Int_t)clov_mult) << " : "
              << No_Clover << " :";
    for (Int_t j = 0; j < No_Clover; j++)
      std::cout << " ( " << ((Int_t)(Clover[j])) << " : " << CLT[j] << " , " << Eclab[j] << " )";
    std::cout << std::endl;
#else /* 0 or 1 */
    // dump "non-zero" entries only
    // if (clov_mult == 0) continue;
    std::cout << i << " : "
              << ((Int_t)clov_mult) << " : "
              << No_Clover << " :";
    for (Int_t j = 0; j < No_Clover; j++)
      if ((CLT[j] != 0) || (Eclab[j] != 0))
        std::cout << " ( " << ((Int_t)(Clover[j])) << " : " << CLT[j] << " , " << Eclab[j] << " )";
    std::cout << std::endl;
#endif /* 0 or 1 */
  }
  
  t->ResetBranchAddresses(); // "disconnect" from local variables
  delete f; // automatically deletes "t", too
  return;
}
