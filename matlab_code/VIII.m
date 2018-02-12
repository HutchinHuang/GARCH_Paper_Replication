% Main Program to Get Table VII.
alphainit = 0.05;
betainit = 0.9;

indices={'ICGN.csv','MCBF.csv','KVA.csv','BTC.csv','CCME.csv'};

[tab_n, tab_alpha, tab_beta, tab_stat, tab_p] = test_the_last(indices,alphainit,betainit);

% disp([tab_n;tab_alpha;tab_beta;tab_stat;tab_p]);
disp(tab_alpha');
disp(tab_beta');
disp(tab_stat');
disp(tab_p');