
package rs.ac.bg.etf.aor2.replacementpolicy;

import rs.ac.bg.etf.aor2.memory.MemoryOperation;
import rs.ac.bg.etf.aor2.memory.cache.ICacheMemory;
import java.util.*;
import rs.ac.bg.etf.aor2.memory.cache.Tag;


public class SLRUReplacementPolicy implements IReplacementPolicy{
    
    protected int protectedSegmentSize;
    protected ICacheMemory cache;
    
    //za svaki set niz zasticenih i mogucih blokova
    protected ArrayList <LinkedList <Integer>> protectedSegments;
    protected ArrayList <LinkedList <Integer>> possibleSegments;
    
    protected int blocksPerSet;
    protected int sets;
    
    @Override
    public void init(ICacheMemory cacheMemory) {
            if(cacheMemory==null)
                 throw new RuntimeException(
                    "Cache missing.");
            
            cache=cacheMemory;
            blocksPerSet=(int)(cache.getBlockNum()/cache.getSetNum());
            protectedSegmentSize=(int)(0.8*blocksPerSet);
            sets=(int)cache.getSetNum();
            protectedSegments=new ArrayList(sets);
            possibleSegments=new ArrayList(sets);
            
            
            for (int i = 0; i < sets; i++) {
            protectedSegments.add(i, new LinkedList<Integer>());
            possibleSegments.add(i, new LinkedList<Integer>());
        }
            
            reset();
    }

    @Override
    public int getBlockIndexToReplace(long adr) {
        int set = (int)cache.extractSet(adr);
        LinkedList <Integer> lista = possibleSegments.get(set);
        int pocetak = (int)set*blocksPerSet;
        ArrayList<Tag> tagMemory = cache.getTags();
        
        //ako ima neki koji je prazan izaberi taj
        for (int i = 0; i < lista.size(); i++) {
            int indeks = lista.get(i);
            if(!tagMemory.get(indeks+pocetak).V)
            {
                lista.remove(i);
                lista.addFirst(indeks);
                return indeks + pocetak;
            }
        }
        
        //iz liste mogucih uzmi najmanje koriscen + azuriraj listu
        int broj = lista.getLast();
        lista.removeLast();
        lista.addFirst(broj);
        return lista.getFirst() + pocetak;
    }

    @Override
    public void doOperation(MemoryOperation operation) {
        MemoryOperation.MemoryOperationType opr = operation.getType();
        ArrayList<Tag> tagMemory = cache.getTags();
        
        if ((opr == MemoryOperation.MemoryOperationType.READ)
                || (opr == MemoryOperation.MemoryOperationType.WRITE)) {

            long adr = operation.getAddress();
            int set = (int) cache.extractSet(adr);
            long tagTag = cache.extractTag(adr);
            LinkedList <Integer> listaMogucih = possibleSegments.get(set);
            LinkedList <Integer> listaZasticenih = protectedSegments.get(set);
            int pocetak = (int)set*blocksPerSet;
            
            //provera da li je u zasticenom delu
            for (int i = 0;i<listaZasticenih.size();i++) {
                int broj = listaZasticenih.get(i);
                Tag tag = tagMemory.get(pocetak + broj);
                if(tag.V && tag.tag==tagTag)
                {
                    listaZasticenih.remove(i);
                    listaZasticenih.addFirst(broj);
                    return;
                }
            }
            
            
            //provera da li je u delu mogucih segmenata
            for (int i = 0;i<listaMogucih.size();i++) {
                int broj = listaMogucih.get(i);
                Tag tag = tagMemory.get(pocetak + broj);
                if(tag.V && tag.tag==tagTag)
                {
                    listaMogucih.remove(i);
                    listaZasticenih.addFirst(broj);
                    if(listaZasticenih.size()>protectedSegmentSize)
                    {
                        int tmp = listaZasticenih.getLast();
                        listaZasticenih.removeLast();
                        listaMogucih.addFirst(tmp);
                    }
                    return;
                }
            }
            
        } else if (operation.getType() == MemoryOperation.MemoryOperationType.FLUSHALL) {
            reset();
            /*for (int i = 0; i < tagMemory.size(); i++) {
                tagMemory.get(i).V = false;
            }*/
        }
           
    }

    @Override
    public String printValid() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String printAll() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void reset() {
            
            for (int i = 0; i < sets; i++) {
            protectedSegments.get(i).clear();
            LinkedList lista = possibleSegments.get(i);
            lista.clear();
                for (int j = 0; j < blocksPerSet; j++) {
                    lista.addLast(j);
                }
        }   
    }
    
}
