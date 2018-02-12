% Main Program to Get Table VII.
alphainit = 0.05;
betainit = 0.9;

indices= {'cac.csv','dax.csv','dja.csv','dowjones.csv','djt.csv', 'dju.csv','ftse.csv','nasdaq.csv','nikkei.csv','smi.csv','sp500.csv'};

res = test_series(indices,alphainit,betainit);
disp(res);


