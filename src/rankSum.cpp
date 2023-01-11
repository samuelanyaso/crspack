// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <cmath>

using namespace std;
using namespace Rcpp;


// [[Rcpp::export]]
bool foo7v2C (arma::mat dw){
	
	arma::vec grps;
	arma::mat dw0;
	arma::mat dw1;
	arma::mat dw0i; 
	arma::mat dw1i;
	arma::uword xcid;
	double nj0; double nj1;
	arma::vec cids = sort(unique(dw.col(0)));
	double m = cids.size();
	arma::vec res(m);
				
	grps = dw.col(2); // group indicators
	dw0 = dw.rows(find(grps == 0)); // subset the data for grp=0
	dw1 = dw.rows(find(grps == 1)); // subset the data for grp=1
	
	for(arma::uword j = 0; j < m; j++) {
		xcid = cids(j);
		
		// for both grps, subset the data and then the outcomes corresponding to cluster i
		dw0i = dw0.rows(find(dw0.col(0) == xcid));
		nj0 = dw0i.n_rows; // ICG size for grp=0

		dw1i = dw1.rows(find(dw1.col(0) == xcid));
		nj1 = dw1i.n_rows; // ICG size for grp=1
		
		if((nj0==0) || (nj1==0)){
			res(j) = 1;
			break;
		}
	}	
	return(any(res == 1));
}

/* Evaluates the third summation in the test statistic - for each xik, returns a vector of length m where the element of j is 0. */
arma::vec gooC (arma::vec dv, arma::mat dw){
	
	arma::vec grps = dw.col(2); // group indicators
	arma::mat dw0 = dw.rows(find(grps == 0)); // subset the data for grp=0
	arma::mat dw1 = dw.rows(find(grps == 1)); // subset the data for grp=1
	
	arma::vec x0 = dw0.col(1); // outcome values for grp=0
	arma::vec x1 = dw1.col(1); // outcome values for grp=0
	arma::vec id0 = dw0.col(0); // cluster ids for patients with grp=0
	arma::vec id1 = dw1.col(0); // cluster ids for patients with grp=1

	arma::vec cids = sort(unique(dw.col(0)));
	double m = cids.size();
	double nj0;
	double nj1;

	double xik = dv(1);
	
	// create a vector to store the results
	arma::vec sl0(m);
	arma::vec sl1(m);
	arma::vec sl2(m);
	
	arma::vec indj = arma::zeros(m); // initialize a zero vector of size m
	
	arma::uword xcid;
	
	arma::vec Fja; arma::vec Fjb; arma::vec Fj;
	
	arma::mat dw0i; arma::vec x0ia; arma::vec x0ib;
	arma::mat dw1i; arma::vec x1ia; arma::vec x1ib;
	arma::mat dwi; arma::vec xia; arma::vec xib;
	
	arma::vec v0; 
	arma::vec v1; 
	arma::vec v2; 
	
	arma::uvec idx0a; arma::uvec idx0b; 
	arma::uvec idx1a; arma::uvec idx1b;
	arma::uvec idx2a; arma::uvec idx2b;
	
	for(arma::uword i = 0; i < m; i++) {
		xcid = cids(i);
		
		// for both grps, subset the data and then the outcomes corresponding to cluster i
		dw0i = dw0.rows(find(dw0.col(0) == xcid));
		x0ia = dw0i.col(1);
		x0ib = x0ia;
		nj0 = x0ia.size(); // ICG size for grp=0
		
		dw1i = dw1.rows(find(dw1.col(0) == xcid));
		x1ia = dw1i.col(1);
		x1ib = x1ia;
		nj1 = x1ia.size(); // ICG size for grp=1
		
		dwi = dw.rows(find(dw.col(0) == xcid));
		xia = dwi.col(1);
		xib = xia;
		
		// for group 0
		if (nj0 != 0) { // the function is written such that the incomplete ICG can only occur in grp0
			idx0a = find(x0ia <= xik);
			x0ia.zeros();
			x0ia.elem(idx0a).fill(1);
			
			idx0b = find(x0ib < xik);
			x0ib.zeros();
			x0ib.elem(idx0b).fill(1);
			v0 = x0ia + x0ib;
				
			sl0(i) = mean(v0) * 0.5;
		} else {
			sl0(i) = 0;
		}
		
		// for group 1
		idx1a = find(x1ia <= xik);
		x1ia.zeros();
		x1ia.elem(idx1a).fill(1);

		idx1b = find(x1ib < xik);
		x1ib.zeros();
		x1ib.elem(idx1b).fill(1);
		v1 = x1ia + x1ib;
		
		sl1(i) = mean(v1) * 0.5;		
		
		if((nj0 > 0) & (nj1 > 0)) {
			indj(i) = 1;
		}	
		
		// for both groups	
		idx2a = find(xia <= xik);
		xia.zeros();
		xia.elem(idx2a).fill(1);
		
		idx2b = find(xib < xik);
		xib.zeros();
		xib.elem(idx2b).fill(1);		
		v2 = xia + xib;
		
		sl2(i) = mean(v2);
	}
	
	Fja = (sl0 + sl1) % indj; // elementwise multiplication of the two vectors
	Fjb = sl2 % (1 - indj); // elementwise multiplication of the two vectors
	Fj = Fja + Fjb;
	
	arma::uvec ik = find(cids == dv(0)); // retrieve the positional index
	Fj.elem(ik).fill(0.0); // replace the value of the vector at the index with 0

	/*
	List res;
	res["ik"] = find(cids == dv(0));
	res["Fj"] = Fj;
	*/
	return(Fj);
}

// [[Rcpp::export]]
/* Evaluates the third summation in the test statistic - for all xijk in data with grp==1.
	Inputs: dw - full data matrix for both grps and all clusters, with three cols - id, outcome, grp.
			dw1 - data matrix for only the grp=1, with two cols - id, outcome.
			m - number of unique clusters
	Output: matrix with nrows of dw1 and ncols is m. */
arma::mat goo2C (arma::mat dw, arma::mat dw1, int m){
	
	// stores the results of the double summation as a matrix where the ncols is the no. of clusters and nrows is the sum of the ICG size for grp 1.
	arma::mat res1(dw1.n_rows, m); 
	arma::vec dv;
	arma::rowvec a3;

	// loop over the cluster units with grp=1
	for (arma::uword k = 0; k < dw1.n_rows; k++){
		
		dv = dw1.row(k).t(); // gets a row vector then convert to a column vector (by default arma::vec is a colvec)
		a3 = gooC(dv, dw).t();
		
		res1.row(k) = a3;
	}
	return(res1);
}


// [[Rcpp::export]]
List DDstatV2C (arma::mat dw, bool iICGs=true, bool get_var=true){
	
	arma::vec cids = sort(unique(dw.col(0)));
	double m = cids.size();
	double nj0;
	double nj1;
	// double ni;
	arma::vec grps;
	arma::mat dw0;
	arma::mat dw1;
	arma::mat dw0i; 
	arma::mat dw1i;
	arma::mat dwi; 
	arma::uword xcid;

	if(iICGs){
		// if the ICG size for group 1 is 0, swap the data corresponding to grp=0 with grp=1.
		// if any n_{i1}=0, that is, such cluster has only patients in grp 0.
		// set the grp indicator for all patients in such cluster to 1.
		// We follow this approach bcos the code is written in terms of n_{i1}.
		// Swapping the grp indicator within the cluster with incomplete ICG size is
		// valid because, resampling a patient from such cluster is not dependent
		// on the grp, since the cluster comprise of only one group.
		
		arma::mat tmp3; 
		grps = dw.col(2); // group indicators
		dw0 = dw.rows(find(grps == 0)); // subset the data for grp=0
		dw1 = dw.rows(find(grps == 1)); // subset the data for grp=1
		
		for(arma::uword j = 0; j < m; j++) {
			xcid = cids(j);
			
			// for both grps, subset the data and then the outcomes corresponding to cluster i
			dw0i = dw0.rows(find(dw0.col(0) == xcid));
			
			dw1i = dw1.rows(find(dw1.col(0) == xcid));
			nj1 = dw1i.n_rows; // ICG size for grp=1
			
			if(nj1 == 0){
				tmp3 = dw1i; 
				
				dw1i = dw0i;
				dw1i.col(2).ones(); // replaces all elements in the third column with 1
				
				dw0i = tmp3;
				
				dwi = join_vert(dw1i,dw0i);
				dw.rows(find(dw.col(0) == xcid)) = dwi; // updates the ith cluster in the full data
			}
		}	
	}
	
	// ****************************************** //
	// Begins computation for the test statistic
	// ****************************************** //
	
	// update the data sets corresponding to both groups
	grps = dw.col(2); // group indicators
	dw0 = dw.rows(find(grps == 0)); // subset the data for grp=0
	dw1 = dw.rows(find(grps == 1)); // subset the data for grp=1
		
	// stores the results of the double summation as a matrix where the ncols is the no. of clusters and nrows is the sum of the ICG size for grp 1.
	double K = dw1.n_rows;
	arma::mat res1(K, cids.size()); 
	arma::vec tmp0(K);
	arma::vec tmp1(K);
	arma::vec tmp2(K);
	arma::vec nv1(K);
	arma::vec tmp4;
	arma::vec tmp5;
	
	double test_stat;
	// for the expected value
	arma::vec res2(cids.size()); 
	double exp_val;
	
	arma::uword cntr = 0; // counter for the nrows in res1
	arma::vec dv;

	// loop over the clusters
	for(arma::uword j = 0; j < m; j++) {
		xcid = cids(j);
		
		// for the entire data set
		dwi = dw.rows(find(dw.col(0) == xcid));
		// ni = dwi.n_rows; // cluster size
		
		// for both grps, subset the data and then the outcomes corresponding to cluster i
		dw0i = dw0.rows(find(dw0.col(0) == xcid));
		nj0 = dw0i.n_rows; // ICG size for grp=0
		
		dw1i = dw1.rows(find(dw1.col(0) == xcid));
		nj1 = dw1i.n_rows; // ICG size for grp=1
		
		// get the first two cols of the data with grp=1
		dw1i = dw1i.cols(0, 1);
		
		// begin computation for the several terms in the formula
		double a1 = 0;
		arma::rowvec a3;
		
		// loop over the cluster units with grp=1
		for (arma::uword k = 0; k < nj1; k++){
			
			if((nj1 > 0) & (nj0 > 0)){
				tmp0(cntr) = 1;
			} else {
				tmp0(cntr) = 0;
			}
			if(nj0 == 0){
				tmp1(cntr) = 1;
			} else {
				tmp1(cntr) = 0;
			};
			nv1(cntr) = nj1;
			
			// for expected value
			a1 = tmp0(cntr) + 2.0 * tmp1(cntr);
			
			dv = dw1i.row(k).t(); // gets a row vector then convert to a column vector (by default arma::vec is a colvec)
			a3 = gooC(dv, dw).t();
			res1.row(cntr) = a3;
			
			cntr++;
		}
		// Computation for the expected value
		// res2(j) = (a1 + 2.0 * a2);
		res2(j) = a1;
	}
	
	exp_val = 0.25*(m + 1)*sum(res2); // expected value of the test statistic
	
	tmp4 = (tmp0 + 2*nv1 % tmp1) / (4 * nv1);
	tmp5 = (tmp0 + 2*tmp1) / (2 * nv1);
		
	// compute the DD statistic
    for (arma::uword k = 0; k < tmp4.size(); k++) {
        res1.row(k) *= tmp4(k);
    }
	test_stat = accu(res1) + accu(tmp5);
	
	// computes the variance
	double vDD; double vDDa; double Zscore;
	arma::vec Tstar;
	if(get_var){
		
		arma::mat res1b;
		arma::vec tmp5b;
		arma::uvec tmp6;
		arma::vec Tstar;
		arma::vec tmp7(m);	
		
		// the leave one cluster out estimates
		for(arma::uword j = 0; j < m; j++) {
			xcid = cids(j);
			tmp6 = find(dw1.col(0) == xcid);
			
			res1b = res1;
			tmp5b = tmp5;
			
			res1b.shed_col(j);
			res1b.shed_rows(tmp6);
			
			tmp5b.shed_rows(tmp6);
			tmp7(j) = accu(res1b) + accu(tmp5b);
		}
		Tstar = test_stat - tmp7;
		vDDa = arma::var(Tstar);
		vDD = ((m*m)/(m-1)) * vDDa;
		if (vDD <= 0.00000001) {
			vDD = 0.00000001;
		} 		
		// Tstar = m*test_stat - (m - 1)*tmp7;
		// vDDa = arma::var(Tstar);
		// vDD = (1/m) * vDDa;
		// if (vDD <= 0.00000001) {
		// 	vDD = 0.00000001;
		// } 	
		Zscore = (test_stat - exp_val)/sqrt(vDD);
	}

	List bes;
	bes["test_stat"] = test_stat;
	bes["exp_val"] = exp_val;
	bes["var"] = vDD;
	bes["Zscore"] = Zscore;
	bes["Tstar"] = Tstar;
	// bes["res1"] = res1;
	// bes["tmp5"] = tmp5;
	// bes["vDDa"] = vDDa;
	return(bes);
}

// [[Rcpp::export]]
List DDstatV1C (arma::mat dw){
	
	arma::vec cids = sort(unique(dw.col(0)));
	double m = cids.size();
	double nj0;
	double nj1;
	double ni;
	arma::mat dw0i; arma::vec x0ia;
	arma::mat dw1i; arma::vec x1ia; 
	arma::mat tmp3; 
	arma::mat dwi; 
	arma::uword xcid;
	arma::mat tmp; 
	
	// if the ICG size for group 1 is 0, swap the data corresponding to grp=0 with grp=1.
	// if any n_{i1}=0, that is, such cluster has only patients in grp 0.
	// set the grp indicator for all patients in such cluster to 1.
	// We follow this approach bcos the code is written in terms of n_{i1}.
	// Swapping the grp indicator within the cluster with incomplete ICG size is
	// valid because, resampling a patient from such cluster is not dependent
	// on the grp, since the cluster comprise of only one group.
	
	arma::vec grps = dw.col(2); // group indicators
	arma::mat dw0 = dw.rows(find(grps == 0)); // subset the data for grp=0
	arma::mat dw1 = dw.rows(find(grps == 1)); // subset the data for grp=1
	
	for(arma::uword j = 0; j < m; j++) {
		xcid = cids(j);
		
		// for both grps, subset the data and then the outcomes corresponding to cluster i
		dw0i = dw0.rows(find(dw0.col(0) == xcid));
		
		dw1i = dw1.rows(find(dw1.col(0) == xcid));
		nj1 = dw1i.n_rows; // ICG size for grp=1
		
		if(nj1 == 0){
			tmp3 = dw1i; 
			
			dw1i = dw0i;
			dw1i.col(2).ones(); // replaces all elements in the third column with 1
			
			dw0i = tmp3;
			
			dwi = join_vert(dw1i,dw0i);
			dw.rows(find(dw.col(0) == xcid)) = dwi; // updates the ith cluster in the full data
		}
	}
	
	// ****************************************** //
	// Begins computation for the test statistic
	// ****************************************** //
	
	// update the data sets corresponding to both groups
	grps = dw.col(2); // group indicators
	dw0 = dw.rows(find(grps == 0)); // subset the data for grp=0
	dw1 = dw.rows(find(grps == 1)); // subset the data for grp=1
	
	// stores the results of the double summation as a matrix where the ncols is the no. of clusters and nrows is the sum of the ICG size for grp 1.
	arma::mat res1(dw1.n_rows, cids.size()); 
	double test_stat;
	// for the expected value
	arma::vec res2(cids.size()); 
	double exp_val;
	
	arma::uword cntr = 0; // counter for the nrows in res1
	arma::vec dv;

	// loop over the clusters
	for(arma::uword j = 0; j < m; j++) {
		xcid = cids(j);
		
		// for the entire data set
		dwi = dw.rows(find(dw.col(0) == xcid));
		ni = dwi.n_rows; // cluster size
		
		// for both grps, subset the data and then the outcomes corresponding to cluster i
		dw0i = dw0.rows(find(dw0.col(0) == xcid));
		nj0 = dw0i.n_rows; // ICG size for grp=0
		
		dw1i = dw1.rows(find(dw1.col(0) == xcid));
		nj1 = dw1i.n_rows; // ICG size for grp=1
		
		// get the first two cols of the data with grp=1
		dw1i = dw1i.cols(0, 1);
		
		// begin computation for the several terms in the formula
		double a1=0; double a2=0; double a3 = 0;
		
		if((nj1 > 0) & (nj0 > 0)){
			a1 = 1;
		} else {
			a1 = 0;
		}
		// if((nj1 > 0) & (nj0 == 0)){
		if(nj0 == 0){
			a2 = 1;
		} else {
			a2 = 0;
		};
		
		tmp = dw;
		tmp.shed_rows(find(dw.col(0) == xcid));
		
		// loop over the cluster units with grp=1
		for (arma::uword k = 0; k < nj1; k++){
			// a3 = 999;
			// Rcpp::Rcout << j << std::endl; // prints out the matrix
			// Rcpp::Rcout << a3 << std::endl; // prints out the matrix
			
			dv = dw1i.row(k).t(); // gets a row vector then convert to a column vector (by default arma::vec is a colvec)
			
			// Rcpp::Rcout << k << std::endl; // prints out the matrix
			// Rcpp::Rcout << dv << std::endl; // prints out the matrix
			
			a3 = sum(gooC(dv, tmp));
			// a3 = sum(gooC(dv, dw));
			
			res1(cntr, j) = (1/(2*nj1)) * (a1 + 2*a2 + ((a1 + 2*(nj1/ni)*a2)/2)*a3);
			cntr++;
		}
		
		// Computation for the expected value
		res2(j) = (a1 + 2*a2);
	}
	
	test_stat = accu(res1); // final test statistic
	exp_val = 0.25*(m + 1)*sum(res2); // expected value of the test statistic
	
	List bes;
	bes["test_stat"] = test_stat;
	bes["exp_val"] = exp_val;

	return(bes);
}


// [[Rcpp::export]]
double var_funC (arma::mat dw){
	
	arma::vec cids = sort(unique(dw.col(0)));
	double m = cids.size();
	arma::vec Tstar(cids.size()); // initialize the vector
	
	List tstat1 = DDstatV1C(dw);
	List tstat3;
	double tstat2 = tstat1["test_stat"];
	double tstat4;
	arma::uword xcid;
	arma::mat tmp;
	double vtest;
	
	for (arma::uword j = 0; j < m; j++) {
		xcid = cids(j);
		tmp = dw;
		tmp.shed_rows(find(dw.col(0) == xcid));
		tstat3 = DDstatV1C(tmp);
		tstat4 = tstat3["test_stat"];
		Tstar(j) = tstat2 - tstat4;
	}
	vtest = ((m*m)/(m-1)) * var(Tstar);
	
	return(vtest);
}


// [[Rcpp::export]]
List ASDstatV1C (arma::mat dw){
	
	// arma::vec Z = sort(dw.col(2)); // the continuous covariate
	arma::vec Z = dw.col(2); // the continuous covariate
	double N = Z.size(); // total number of observations
	double z;
	arma::mat tmp1; // new matrix where third col. will contain the group labels.
	arma::vec tmp3;
	arma::uvec indx;
	List tstat1; // stores the results for the T stat
	double tmp2; // stores the result for the squared difference.
	arma::vec Vstar(Z.size()); // initialize the vector
	double v_stat; double a1; double a2;
	
	for (arma::uword j = 0; j < N; j++) {
		// Rcpp::Rcout << j << std::endl; // prints out iteration number

		z = Z(j);
		tmp1 = dw;
		tmp3 = tmp1.col(2);
		tmp3.zeros();
		indx = find(dw.col(2) <= z);
		tmp3.elem(indx).ones();
		tmp1.col(2) = tmp3;
				
		tstat1 = DDstatV1C(tmp1);
		a1 = tstat1["test_stat"]; a2 = tstat1["exp_val"];
		tmp2 = a1 - a2;
		Vstar(j) = tmp2 * tmp2;
	}
	v_stat = accu(Vstar)/N;
	
	List res;
	res["v_stat"] = v_stat;
	res["Vstar"] = Vstar;
	// return(v_stat);
	return(res);
}

// [[Rcpp::export]]
List ASDstatV2C (arma::mat dw, bool iICGs){
  
	// arma::vec Z = sort(dw.col(2)); // the continuous covariate
	arma::vec Z = dw.col(2); // the continuous covariate
	double N = Z.size(); // total number of observations
	double z;
	arma::mat tmp1; // new matrix where third col. will contain the group labels.
	arma::vec tmp3;
	arma::uvec indx;
	List tstat1; // stores the results for the T stat
	double tmp2; // stores the result for the squared difference.
	arma::vec Vstar(Z.size()); // initialize the vector
	double v_stat; 
	// double a1; double a2; double a3;
	bool b1; bool b2;
	double N0 = 0;
	
	for (arma::uword j = 0; j < N; j++) {
		// Rcpp::Rcout << j << std::endl; // prints out iteration number

		z = Z(j);
		tmp1 = dw;
		tmp3 = tmp1.col(2);
		tmp3.zeros();
		indx = find(Z <= z);
		tmp3.elem(indx).ones();
		tmp1.col(2) = tmp3;
		
		if(iICGs){
			b1 = foo7v2C(tmp1);
			b2 = false;
		} else {
			b1 = false;
			b2 = true;
		}
		
		// move to the next z if dataset with incomplete ICG str is created
		// else, compute the DD statistic
		if(b1){
			N0++;
			continue;
		} else {
			tstat1 = DDstatV2C(tmp1, b2, true);
			// a1 = tstat1["test_stat"]; a2 = tstat1["exp_val"]; a3 = tstat1["var"];
			// tmp2 = (a1 - a2)/sqrt(a3);
			tmp2 = tstat1["Zscore"]; 
			Vstar(j) = tmp2 * tmp2;
		}
	}
	v_stat = accu(Vstar)/(N - N0);
	
	List res;
	res["v_stat"] = v_stat;
	res["Vstar"] = Vstar;
	return(res);
}


// [[Rcpp::export]]
List ASDpvC (arma::mat dw, int K){
	
	arma::vec cids = sort(unique(dw.col(0))); // unique cluster ids
	double m = cids.size(); // total number of unique clusters
	arma::uword xcid;
	List tmp;
	arma::mat tmp1;
	arma::vec tmp2;
	arma::vec tmp3;
	NumericVector tmp4a; // Rcpp NumericVector
	arma::vec tmp4b; 
	arma::uvec indx;
	arma::vec vres(K); // initialize the vector
	arma::vec tmp5;
	double p_value;
	double obsV; 
	
	// obtain the observed test statistic
	tmp = ASDstatV1C(dw); 
	obsV = tmp["v_stat"]; 
	
	for (int k = 0; k < K; k++) {
		tmp1 = dw;
		tmp2 = tmp1.col(2); // continuous covariate for all clusters
		
		// permute the Zs within each cluster
		for (arma::uword j = 0; j < m; j++) {
			xcid = cids(j);
			indx = find(dw.col(0) == xcid);
			tmp3 = tmp2.elem(indx); // covariate corresponding to cluster j
			tmp4a = RcppArmadillo::sample(tmp3, tmp3.size(), false); // similar to R's sample function
			tmp4b = as<arma::vec>(tmp4a);   // convert Rcpp::NumericVector to arma::vec
			tmp2.elem(indx) = tmp4b; // replaces the original values of Z in cluster j with the permuted values
		}
		tmp1.col(2) = tmp2; // replace the Z column with permuted vector
		
		// compute the V statistic for the permuted sample
		tmp = ASDstatV1C(tmp1);
		vres(k) = tmp["v_stat"]; 
	}
	
	// calculate the p-value
	tmp5.zeros(vres.size());
	tmp5.elem(find(vres >= obsV)).fill(1.0);
	p_value = accu(tmp5)/K;
	
	List res;
	res["obsV"] = obsV;
	res["vres"] = vres;
	res["p_value"] = p_value;
	
	return(res);
}

arma::vec foo6v2 (arma::mat dw){
	
	// arma::vec Z = sort(dw.col(2)); // the continuous covariate
	arma::vec Z = dw.col(2); // the continuous covariate
	double N = Z.size(); // total number of observations
	double z;
	arma::vec grps;
	arma::mat dw0;
	arma::mat dw1;
	arma::mat dw0i; 
	arma::mat dw1i;
	arma::uword xcid;
	double nj0; double nj1;
	arma::mat tmp1; // new matrix where third col. will contain the group labels.
	arma::vec tmp3;
	arma::uvec indx;
	arma::vec cids = sort(unique(dw.col(0)));
	double m = cids.size();
	arma::vec res(N);
	// std::vector<bool> res2(N); // initializes a boolean vector of 
	
	for (arma::uword k = 0; k < N; k++) {
		z = Z(k);
		tmp1 = dw;
		tmp3 = tmp1.col(2);
		tmp3.zeros();
		indx = find(dw.col(2) <= z);
		tmp3.elem(indx).ones();
		tmp1.col(2) = tmp3;
				
		grps = tmp1.col(2); // group indicators
		dw0 = tmp1.rows(find(grps == 0)); // subset the data for grp=0
		dw1 = tmp1.rows(find(grps == 1)); // subset the data for grp=1
		
		for(arma::uword j = 0; j < m; j++) {
			xcid = cids(j);
			
			// for both grps, subset the data and then the outcomes corresponding to cluster i
			dw0i = dw0.rows(find(dw0.col(0) == xcid));
			nj0 = dw0i.n_rows; // ICG size for grp=0

			dw1i = dw1.rows(find(dw1.col(0) == xcid));
			nj1 = dw1i.n_rows; // ICG size for grp=1
			
			if((nj0==0) || (nj1==0)){
				res(k) = 1;
				break;
			}
		}	
	}
	return(res);
}