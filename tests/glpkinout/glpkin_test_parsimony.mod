/* sets */
set PROTEINS;


/* decision variables: yi, i in {1,..,5}. yi = 1 -> protein i is selected */
var y {i in PROTEINS} binary >=0;
/* objective function */
minimize z: sum{i in PROTEINS} y[i];

/* Constraints */
s.t. c1: y[6] >=1;
s.t. c2: y[3] >=1;
s.t. c3: y[2] >=1;
s.t. c4: y[3] >=1;
s.t. c5: y[3] >=1;
s.t. c6: y[2] >=1;
s.t. c7: y[0]+y[1] >=1;
s.t. c8: y[3] >=1;
s.t. c9: y[3] >=1;
s.t. c10: y[3] >=1;
s.t. c11: y[6]+y[7] >=1;
s.t. c12: y[4] >=1;
s.t. c13: y[5] >=1;
s.t. c14: y[3] >=1;
s.t. c15: y[3] >=1;
s.t. c16: y[0]+y[1]+y[2] >=1;
s.t. c17: y[3] >=1;
s.t. c18: y[0]+y[2] >=1;
s.t. c19: y[0] >=1;
s.t. c20: y[3] >=1;
s.t. c21: y[1]+y[2] >=1;
s.t. c22: y[3] >=1;
s.t. c23: y[3] >=1;
s.t. c24: y[2] >=1;
s.t. c25: y[4] >=1;

data;
set PROTEINS := 5 0 6 1 7 2 3 4 ;

end;