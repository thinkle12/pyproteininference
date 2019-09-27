/* sets */
set PROTEINS;


/* decision variables: yi, i in {1,..,5}. yi = 1 -> protein i is selected */
var y {i in PROTEINS} binary >=0;
/* objective function */
minimize z: sum{i in PROTEINS} y[i];

/* Constraints */
s.t. c1: y[7] >=1;
s.t. c2: y[2]+y[6] >=1;
s.t. c3: y[3] >=1;
s.t. c4: y[3] >=1;
s.t. c5: y[3] >=1;
s.t. c6: y[1]+y[7]+y[3] >=1;
s.t. c7: y[4] >=1;
s.t. c8: y[0] >=1;
s.t. c9: y[7]+y[3] >=1;
s.t. c10: y[2] >=1;
s.t. c11: y[1]+y[3] >=1;
s.t. c12: y[0]+y[5] >=1;
s.t. c13: y[1] >=1;
s.t. c14: y[1]+y[7] >=1;

data;
set PROTEINS := 4 1 7 0 5 3 2 6 ;

end;