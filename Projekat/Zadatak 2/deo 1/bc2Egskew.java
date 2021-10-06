package rs.ac.bg.etf.predictor.bc2egskew;

import rs.ac.bg.etf.automaton.Automaton;
import rs.ac.bg.etf.predictor.BHR;
import rs.ac.bg.etf.predictor.Instruction;
import rs.ac.bg.etf.predictor.Predictor;
import rs.ac.bg.etf.predictor.bimodal.Bimodal;
import rs.ac.bg.etf.predictor.egskew.EGskew;

public class bc2Egskew implements Predictor {

    Automaton[] bank1, bank2, meta;
    Bimodal bank0;
    BHR bhr;
    int updateNum;
    boolean bimodalPred, metaPred, egskewPred;

    public bc2Egskew(int BimodalBHRsize, int BimodalnumOfLastBitsAddressGroups, int BimodalnumOfLastBitsAddressSelector, Automaton.AutomatonType Bimodaltype,
            int BHRSize, Automaton.AutomatonType type, int phtSize) {

        bank0 = new Bimodal(BimodalBHRsize, BimodalnumOfLastBitsAddressGroups, BimodalnumOfLastBitsAddressSelector, Bimodaltype);
        bhr = new BHR(BHRSize);
        bank1 = Automaton.instanceArray(type, 1 << phtSize);
        bank2 = Automaton.instanceArray(type, 1 << phtSize);
        meta = Automaton.instanceArray(type, 1 << phtSize);
        updateNum = 0;
        metaPred = bimodalPred = egskewPred = false;

    }

    @Override
    public boolean predict(Instruction branch) {
        boolean metaPredict = meta[hash3(branch)].predict();

        boolean bimodalPredict = bank0.predict(branch);
        boolean bank1Predict = bank1[hash1(branch)].predict();
        boolean bank2Predict = bank2[hash2(branch)].predict();
        if (!(bank1Predict ^ bank2Predict^bimodalPredict)) {
            updateNum = 3;
            egskewPred = bank2Predict;
        } else if (!(bimodalPredict ^ bank2Predict)) {
            updateNum = 1;
            egskewPred = bimodalPredict;
        } else if (!(bimodalPredict ^ bank1Predict)) {
            updateNum = 0;
            egskewPred = bank1Predict;
        } else {
            updateNum = 2;
            egskewPred = bank2Predict;
        }
        if (metaPredict) {
            metaPred = true;
            return bank0.predict(branch);
        } else {
            metaPred = false;
            return egskewPred;
        }
    }

    @Override
    public void update(Instruction branch) {
        boolean outcome = branch.isTaken();
        bhr.insertOutcome(outcome);
        if (metaPred) {
            bank0.update(branch);
        } else {
            switch (updateNum) {
                case 0:
                    bank0.update(branch);
                    bank1[hash1(branch)].updateAutomaton(outcome);
                    break;
                case 1:
                    bank0.update(branch);
                    bank2[hash2(branch)].updateAutomaton(outcome);
                    break;
                case 2:
                    bank2[hash2(branch)].updateAutomaton(outcome);
                    bank1[hash1(branch)].updateAutomaton(outcome);
                    break;
                case 3:
                    bank0.update(branch);
                    bank2[hash2(branch)].updateAutomaton(outcome);
                    bank1[hash1(branch)].updateAutomaton(outcome);
                    break;
                default:
                    break;
            }
        }
        
        if(!(bimodalPred^egskewPred))
        {
            meta[hash3(branch)].updateAutomaton(outcome);
        }

    }

    public Bimodal dohvatiBimodal() {
        return bank0;
    }

    private int hash1(Instruction branch) {
        return ((int) branch.getAddress() ^ bhr.getValue()) % bank1.length;
    }

    private int hash2(Instruction branch) {
        return ((int) branch.getAddress() & (bhr.getValue() >> 2)) % bank2.length;
    }

    private int hash3(Instruction branch) {
        return ((int) branch.getAddress() ^ (bhr.getValue() >> 2)) % bank2.length;
    }
}
