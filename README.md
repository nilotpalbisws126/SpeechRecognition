# Digit Recognition System

**Author:** Nilotpal Biswas  
**Roll No.:** 176101101  
**Date:** 20th Nov 2018

---

## Getting Started

1. Clone or download this repository.  
2. Open a terminal and navigate to the **`DigitRecognitionAssignment/`** folder.  
3. Build and run the main program:

## **************** File info****************

feature_ext.h has 3 main funtions train(), buildHMM() and testing()
Where train() reads the training digit data and writes a file called ciFeatures.txt after calculating capstrel coefficients.
buildHMM() takes the Ci values from codebook and makes observation sequence and also writes observationSeq.txt file under write folder and also it 
writes the filnal HMM model for each digit under /write/HMM_training folder
testing() takes the voice input and predicts what is said calculating the probability for each digit.

make_codebook.h has lbg() function which makes a codebook from the Ci values and writes a file LBG_log.txt under write folder.


```bash
cd DigitRecognitionAssignment
./DigitRecognitionAssignment    # or run the executable produced by your build system


