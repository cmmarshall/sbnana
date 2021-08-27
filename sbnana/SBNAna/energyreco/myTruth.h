#ifndef myTruth_h
#define myTruth_h

//==== Number of particles
const Var varTruthNNeutron = SIMPLEVAR(truth.nneutron);
const Var varTruthNPiMinus = SIMPLEVAR(truth.npiminus);
const Var varTruthNPiPlus = SIMPLEVAR(truth.npiplus);
const Var varTruthNChargedPion = varTruthNPiMinus+varTruthNPiPlus;
const Var varTruthNPiZero = SIMPLEVAR(truth.npizero);
const Var varTruthNProton = SIMPLEVAR(truth.nproton);

#endif
