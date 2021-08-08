
# BCS Sorting.....

to compile: 
```
cmake3 CMakeList.txt
make 
```

## BCSint.cxx
**DoSort()**: Organize raw data(run.file) to correalted implant and decay(beta.root) 
            Input: run*.file => tree:dchan;
            Output: beta*.root => tree: beta;
                    event*.root => tree: event; 
                    output*.root => histograms including pid and so on;

**TOFfluctuation()**: Get TH2D of TOF as runtime;
                    Input: event*.root;
                    Output: tof*.root including (TH2D) tof*;

**CorrectTOF()**: Correct TOF fluctuation;
                  *Request: tof*.root must exist;
                  Input: event*.root;
                  Output: ctof*.root including (Th2D) ctof*;

