
package rs.ac.bg.etf.predictor.egskew;

import rs.ac.bg.etf.automaton.Automaton;
import rs.ac.bg.etf.predictor.BHR;
import rs.ac.bg.etf.predictor.Instruction;
import rs.ac.bg.etf.predictor.Predictor;
import rs.ac.bg.etf.predictor.bimodal.Bimodal;


public class EGskew implements Predictor {
    
    Automaton[]bank1,bank2;
    Bimodal bank0;
    BHR bhr;
    int updateNum;
    
    public EGskew(int BimodalBHRsize, int BimodalnumOfLastBitsAddressGroups, int BimodalnumOfLastBitsAddressSelector, Automaton.AutomatonType Bimodaltype, 
            int BHRSize, Automaton.AutomatonType type, int phtSize){
        
        bank0= new Bimodal(BimodalBHRsize,BimodalnumOfLastBitsAddressGroups,BimodalnumOfLastBitsAddressSelector,Bimodaltype);
        bhr = new BHR(BHRSize);
        bank1 = Automaton.instanceArray(type, 1 << phtSize);
        bank2 = Automaton.instanceArray(type, 1 << phtSize);
        updateNum=0;
        
    }
    
    @Override
    public boolean predict(Instruction branch) {
        boolean bimodalPredict = bank0.predict(branch);
        boolean bank1Predict = bank1[hash1(branch)].predict();
        boolean bank2Predict = bank2[hash2(branch)].predict();
        if(!(bank1Predict^bank2Predict^bimodalPredict)){
            updateNum = 3;
            return bank2Predict;
        }
         else if(!(bimodalPredict^bank2Predict))
         {
             updateNum = 1;
             return bimodalPredict;
         }
         else if(!(bimodalPredict^bank1Predict)){
            updateNum = 0;
            return bank1Predict;
         }
         else{
             updateNum = 2;
             return bank2Predict;
         }
    }
    
    @Override
    public void update(Instruction branch) {
        boolean outcome = branch.isTaken();
        bhr.insertOutcome(outcome);
        switch(updateNum){
            case 0: bank0.update(branch);
                    bank1[hash1(branch)].updateAutomaton(outcome);
                    break;
            case 1:bank0.update(branch);
                   bank2[hash2(branch)].updateAutomaton(outcome);
                   break;
            case 2:bank2[hash2(branch)].updateAutomaton(outcome);
                   bank1[hash1(branch)].updateAutomaton(outcome);
                   break;
            case 3:bank0.update(branch);
                   bank2[hash2(branch)].updateAutomaton(outcome);
                   bank1[hash1(branch)].updateAutomaton(outcome);
                   break;
            default:break;
        }
       
    }
    
    public Bimodal dohvatiBimodal(){
            return bank0;
    }
    
    private int hash1(Instruction branch)
    {
        return ((int)branch.getAddress()^bhr.getValue())%bank1.length;
    }
    
    private int hash2(Instruction branch)
    {
        return ((int)branch.getAddress()^(bhr.getValue()>>2))%bank2.length;
    }
}
