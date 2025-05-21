import ROOT

ROOT.gInterpreter.Declare('#include "../interface/JetInfoDS.h"')

inputFile = 'TCPNtuple.root'


jets = ROOT.JetInfoDS()

fchain = ROOT.TChain('tcpNtuples/analysisTree')
fchain.Add(inputFile)

fchain.SetBranchAddress("Jets", ROOT.AddressOf(jets))

print(fchain.GetEntries())

for iev in range(fchain.GetEntries()): # Be careful!!!
    fchain.GetEntry(iev)
    print(iev, jets.at(0).pt)
