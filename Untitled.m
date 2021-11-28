for i=11:100
   close all
   load("C:\Users\alber\Desktop\muVES retraining\data2\xtrain\"+string(i)+".mat"); 
   load("C:\Users\alber\Desktop\muVES retraining\data2\ytrain\"+string(i)+".mat");
   figure();volshow(x);figure();volshow(y)
end