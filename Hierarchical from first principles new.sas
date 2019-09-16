options ls=130 ;
data a ;
input X1	X2	X3 ;
cards;
25.431565	25.933835	37.24979
23.337532	33.238455	39.632566
12.946236	31.090874	41.410102
11.969851	20.15965	35.350383
14.848399	22.939167	36.585108
18.096521	22.285654	40.311387
15.983284	35.065821	46.766173
20.580939	28.317143	41.124029
19.602993	26.171693	35.882251
4.3317016	29.345467	37.910276
19.086121	26.295734	30.270269
13.548145	25.079301	37.958955
118.21118	135.7007	144.21129
116.95557	132.1643	139.54306
119.87893	126.80532	140.47414
126.91946	132.77775	139.48535
125.7058	124.3089	134.14785
117.85024	128.38386	136.07865
129.14882	131.11702	141.80763
124.01972	134.28581	138.9426
119.44251	130.94937	139.70347
116.28171	130.44813	146.72788
120.91516	129.20598	142.62789
117.39883	125.84577	133.06846
66.004697	76.209273	86.272148
67.442937	79.967425	87.070665
62.346243	74.295596	85.055967
73.360764	81.442888	90.72273
77.058177	82.016948	87.103163
64.71757	88.460159	93.582392
70.499663	85.34535	90.899389
61.26967	80.649116	92.365147
80.102074	83.83776	91.764146
77.419293	81.371541	90.413738
70.364107	77.508103	87.324323
66.789004	88.558624	92.471582
; 

run ;
* From observation 9 up until the end of row 15;
data a ; 
set a  (firstobs=9 obs=15) ;
run ;

proc print data=a ;
run ;

proc gplot data=a ;
plot  x2*x1 ;
run ;

proc iml ;

use  a ;
read all into x ;
print x ;

n=nrow(x) ;
p=ncol(x) ;
print x n p ;
nclus=2 ;


cols=1:n ;
print cols ;
* character matrix and compress function;
nm= rowcatc(compress(    J(n,1,"o" ) || char(((1:n)`))  ));
nm= nm // "-------------------------------------------------------------------------------------------------------------------------------------------------" ;
print nm ;

* Take only the first 7 rows of nm;
clres = nm[1:nrow(nm)-1] ;
print clres ;

dist=distance(x) ;
mmdist=round(max(dist)+9900) ;
dist= diag(J(n,1,mmdist))+dist ; ;
print dist ;

do i = 1 to (n-nclus); *n=7, nclus=2 so iterate 5 times;
 print "***********************" ;
 print "Iteration: " i ;

* Get the minimum value of each column;
	mind = min(dist) ;
	cmind = dist[><,] ;
* Locate where minimum row vector equals the minimum found for the whole distance matrix;
	posmin = loc(cmind=mind) ;
	print mind cmind posmin;

* Single linkage - Concatenate column 1 and 3 then give the minimum of each row;
	distcom = 
	 (dist[,posmin[1]] || dist[,posmin[2]])[,><] ;
	print distcom ;
* Change the minimum in column and rows of posmin to that max;
	distcom[posmin[1],1] = mmdist ;
	distcom[posmin[2],1] = mmdist ;

	print distcom ;

	dist[posmin[1],] = distcom` ;
	dist[,posmin[1]] = distcom  ;
	
* Locate cols that do not include pomin[2]/column 3;
	ncols = loc(cols^=posmin[2]) ;
	print cols ncols ;
	names =  nm[ncols] ; * In the nm matrix only the columns in ncols(without posmin[2]);
	print names ;
	names[posmin[1]] = rowcatc(nm[posmin[1]] || "+" || nm[posmin[2]]) ; * Where nm is = to posmin[1]/1st row;
	print names ;

* Print everything;
	dist=dist[ncols,ncols] ;
	nm=names ;
	nm= nm // "-------------------------------------------------------------------------------------------------------------------------------------------------" ;
	cols=   1:(n-i) ; * original cols=1:n, ncols ; * As the number of clusters increase(i) the columns of the matrix decrease;
	print "AAA" cols;

	print dist ,  nm;
	
	clres = clres || ( (nm[1:nrow(nm)-1]) // J(i,1,"-0-") );
	print clres  ;
end ;
*end of loop;

cls =   rowcatc( J(n,1,"Cluster= ") || char(n:1)`) ;
print cls ;

print clres[colname=cls] ;

fres = clres[1:nclus,ncol(clres)] ;
print fres ;

quit ;

