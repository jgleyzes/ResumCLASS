#include "ResumCLASS.h"

double I00 (const double & k, const double & X1q, const double & f1) {

	return 1.*pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (2. + 0.022222222222222227*pow(f1,2)*pow(2. + f1,2)*pow(k,4)*pow(X1q,2) - 
       0.000705467372134039*pow(f1,3)*pow(2. + f1,3)*pow(k,6)*pow(X1q,3)) ;
}

double I02 (const double & k, const double & X1q, const double & f1) {

	return f1*X1q*pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*pow(k,2)*
     (-0.26666666666666666 + X1q*pow(f1,3)*pow(k,2)*(0.0063492063492063475 - 0.012698412698412709*X1q*pow(k,2)) + 
       X1q*pow(f1,2)*pow(k,2)*(0.02539682539682539 - 0.008465608465608468*X1q*pow(k,2)) + 
       f1*(-0.13333333333333333 + 0.02539682539682539*X1q*pow(k,2)) - 0.006349206349206354*pow(f1,4)*pow(k,4)*pow(X1q,2) - 
       0.0010582010582010585*pow(f1,5)*pow(k,4)*pow(X1q,2)) ;
}

double I04 (const double & k, const double & X1q, const double & f1) {

	return f1*X1q*pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*pow(k,2)*
     (-2.220446049250313e-16 + X1q*pow(f1,3)*pow(k,2)*(0.006349206349206375 - 0.004617604617604604*X1q*pow(k,2)) + 
       X1q*pow(f1,2)*pow(k,2)*(0.0253968253968255 - 0.0030784030784030973*X1q*pow(k,2)) + 
       f1*(-1.1102230246251565e-16 + 0.0253968253968255*X1q*pow(k,2)) - 0.0023088023088023157*pow(f1,4)*pow(k,4)*pow(X1q,2) - 
       0.00038480038480038364*pow(f1,5)*pow(k,4)*pow(X1q,2)) ;
}

double I06 (const double & k, const double & X1q, const double & f1) {

	return f1*X1q*pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*pow(k,2)*
     (-1.7763568394002503e-15 + X1q*pow(f1,3)*pow(k,2)*(-1.6653345369377348e-16 - 0.002664002664002574*X1q*pow(k,2)) + 
       X1q*pow(f1,2)*pow(k,2)*(-4.4408920985006257e-16 - 0.0017760017760017899*X1q*pow(k,2)) + 
       f1*(-8.881784197001251e-16 - 6.661338147750939e-16*X1q*pow(k,2)) - 0.001332001332001287*pow(f1,4)*pow(k,4)*pow(X1q,2) - 
       0.00022200022200022373*pow(f1,5)*pow(k,4)*pow(X1q,2)) ;
}

double I08 (const double & k, const double & X1q, const double & f1) {

	return f1*X1q*pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*pow(k,2)*
     (-3.552713678800501e-15 - 1.7763568394002505e-15*f1 + 
       X1q*pow(f1,2)*pow(k,2)*(-8.881784197001252e-16 - 2.220446049250313e-16*X1q*pow(k,2)) - 
       1.3322676295501878e-15*pow(f1,3)*pow(k,4)*pow(X1q,2) - 6.661338147750939e-16*pow(f1,4)*pow(k,4)*pow(X1q,2) - 
       1.1102230246251565e-16*pow(f1,5)*pow(k,4)*pow(X1q,2)) ;
}

double I22 (const double & k, const double & X1q, const double & f1) {

	return  pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (0.4000000000000002 - 0.07619047619047603*f1*X1q*pow(k,2) + 
       X1q*pow(f1,2)*pow(k,2)*(-0.038095238095238015 + 0.038095238095238126*X1q*pow(k,2)) + 
       pow(f1,4)*(0.009523809523809532 - 0.007696007696007687*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.038095238095238126 - 0.0051306717973384555*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.0038480038480038434*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.0006413339746673069*pow(f1,6)*pow(k,6)*pow(X1q,3)) ;
}

double I24 (const double & k, const double & X1q, const double & f1) {

	return  pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (-2.220446049250313e-16 - 0.07619047619047647*f1*X1q*pow(k,2) + 
       X1q*pow(f1,2)*pow(k,2)*(-0.03809523809523824 + 0.013852813852813894*X1q*pow(k,2)) + 
       pow(f1,4)*(0.0034632034632034805 - 0.006038406038406163*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.013852813852813922 - 0.00402560402560409*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.0030192030192030817*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.0005032005032005113*pow(f1,6)*pow(k,6)*pow(X1q,3)) ;
}

double I26 (const double & k, const double & X1q, const double & f1) {

	return  pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (1.7763568394002505e-15 + 0.007992007992009054*pow(f1,2)*pow(k,4)*pow(X1q,2) + 
       pow(f1,4)*(0.0019980019980022634 - 0.0021312021312017704*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.007992007992009054 - 0.0014208014208013653*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.0010656010656008852*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.00017760017760017066*pow(f1,6)*pow(k,6)*pow(X1q,3)) ;
}

double I28 (const double & k, const double & X1q, const double & f1) {

	return pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*pow(f1,2)*
     (1.7763568394002505e-15 - 0.0004387769093637528*X1q*pow(f1,3)*pow(k,2) - 0.00007312948489415771*X1q*pow(f1,4)*pow(k,2) + 
       pow(f1,2)*(4.440892098500626e-16 - 0.0008775538187275056*X1q*pow(k,2)) + 
       f1*(1.7763568394002505e-15 - 0.0005850358791532617*X1q*pow(k,2)))*pow(k,4)*pow(X1q,2) ;
}

double I44 (const double & k, const double & X1q, const double & f1) {

	return  pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (0.2222222222222241 - 0.03848003848003678*f1*X1q*pow(k,2) + 
       X1q*pow(f1,2)*pow(k,2)*(-0.01924001924001839 + 0.017651484318150557*X1q*pow(k,2)) + 
       pow(f1,4)*(0.004412871079537639 - 0.004058492947381698*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.017651484318150557 - 0.0027056619649215574*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.002029246473690849*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.0003382077456151947*pow(f1,6)*pow(k,6)*pow(X1q,3)) ;
}

double I46 (const double & k, const double & X1q, const double & f1) {

	return pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (3.552713678800501e-15 - 0.04662004662005259*f1*X1q*pow(k,2) + 
       X1q*pow(f1,2)*pow(k,2)*(-0.023310023310026295 + 0.007992007992005279*X1q*pow(k,2)) + 
       pow(f1,4)*(0.0019980019980013197 - 0.0032647091470612293*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.007992007992005279 - 0.002176472764708448*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.0016323545735306146*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.000272059095588556*pow(f1,6)*pow(k,6)*pow(X1q,3)) ;
}

double I48 (const double & k, const double & X1q, const double & f1) {

	return  pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*pow(f1,2)*
     (0.0051190639425708895 - 0.0006543164437946558*X1q*pow(f1,3)*pow(k,2) - 0.00010905274063233161*X1q*pow(f1,4)*pow(k,2) + 
       pow(f1,2)*(0.0012797659856427224 - 0.0013086328875893116*X1q*pow(k,2)) + 
       f1*(0.0051190639425708895 - 0.0008724219250586529*X1q*pow(k,2)))*pow(k,4)*pow(X1q,2) ;
}

double I66 (const double & k, const double & X1q, const double & f1) {

	return    pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (0.15384615384617462 - 0.026107226107185966*f1*X1q*pow(k,2) + 
       X1q*pow(f1,2)*pow(k,2)*(-0.013053613053592983 + 0.011956670780190137*X1q*pow(k,2)) + 
       pow(f1,4)*(0.0029891676950475343 - 0.0026729834779315453*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.011956670780190137 - 0.0017819889852868087*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.0013364917389657727*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.00022274862316107313*pow(f1,6)*pow(k,6)*pow(X1q,3)) ;
}

double I68 (const double & k, const double & X1q, const double & f1) {

	return  f1*X1q*pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*pow(k,2)*
     (-0.03378582202131497 + X1q*pow(f1,3)*pow(k,2)*(0.0014279494155626793 - 0.0023247478356616114*X1q*pow(k,2)) + 
       X1q*pow(f1,2)*pow(k,2)*(0.005711797662250717 - 0.0015498318904256791*X1q*pow(k,2)) + 
       f1*(-0.016892911010657485 + 0.005711797662250717*X1q*pow(k,2)) - 0.0011623739178308057*pow(f1,4)*pow(k,4)*pow(X1q,2) - 
       0.0001937289863027658*pow(f1,5)*pow(k,4)*pow(X1q,2)) ;
}

double I88 (const double & k, const double & X1q, const double & f1) {

	return pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (0.11764705882364977 - 0.019814241485391904*f1*X1q*pow(k,2) + 
       X1q*pow(f1,2)*pow(k,2)*(-0.009907120742695952 + 0.009075479044724943*X1q*pow(k,2)) + 
       pow(f1,4)*(0.002268869761181236 - 0.002012969357252814*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.009075479044724943 - 0.0013419795712650284*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.001006484678626407*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.00016774744640812855*pow(f1,6)*pow(k,6)*pow(X1q,3)) ;
}

//////

double J00 (const double & k, const double & X1q, const double & f1) {

	return pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (0.6666666666666666 - 0.17777777777777784*f1*X1q*pow(k,2) + 
       X1q*pow(f1,2)*pow(k,2)*(-0.08888888888888892 + 0.04656084656084655*X1q*pow(k,2)) + 
       pow(f1,4)*(0.011640211640211638 - 0.011287477954144622*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.04656084656084655 - 0.007524985302763083*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.005643738977072311*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.0009406231628453854*pow(f1,6)*pow(k,6)*pow(X1q,3)) ;
}

double J02 (const double & k, const double & X1q, const double & f1) {

	return pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (0.26666666666666666 - 0.13968253968253963*f1*X1q*pow(k,2) + 
       X1q*pow(f1,2)*pow(k,2)*(-0.06984126984126982 + 0.03386243386243387*X1q*pow(k,2)) + 
       pow(f1,4)*(0.008465608465608468 - 0.009363476030142693*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.03386243386243387 - 0.006242317353428459*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.0046817380150713465*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.0007802896691785574*pow(f1,6)*pow(k,6)*pow(X1q,3)) ;
}

double J04 (const double & k, const double & X1q, const double & f1) {

	return  f1*X1q*pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*pow(k,2)*
     (-0.050793650793650835 + X1q*pow(f1,3)*pow(k,2)*(0.004425204425204443 - 0.005564805564805625*X1q*pow(k,2)) + 
       X1q*pow(f1,2)*pow(k,2)*(0.017700817700817773 - 0.0037098703765370694*X1q*pow(k,2)) + 
       f1*(-0.025396825396825418 + 0.017700817700817773*X1q*pow(k,2)) - 0.0027824027824028125*pow(f1,4)*pow(k,4)*pow(X1q,2) - 
       0.0004637337970671337*pow(f1,5)*pow(k,4)*pow(X1q,2)) ;
}

double J22 (const double & k, const double & X1q, const double & f1) {

	return pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (0.20952380952380945 - 0.10158730158730161*f1*X1q*pow(k,2) + 
       X1q*pow(f1,2)*pow(k,2)*(-0.05079365079365081 + 0.028090428090428055*X1q*pow(k,2)) + 
       pow(f1,4)*(0.007022607022607014 - 0.007794674461341154*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.028090428090428055 - 0.00519644964089408*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.003897337230670577*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.00064955620511176*pow(f1,6)*pow(k,6)*pow(X1q,3));
}

double J24 (const double & k, const double & X1q, const double & f1) {

	return  pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (0.07619047619047647 - 0.05310245310245332*f1*X1q*pow(k,2) + 
       X1q*pow(f1,2)*pow(k,2)*(-0.02655122655122666 + 0.016694416694416736*X1q*pow(k,2)) + 
       pow(f1,4)*(0.004173604173604184 - 0.005170138503471949*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.016694416694416736 - 0.0034467590023146144*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.0025850692517359747*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.0004308448752893268*pow(f1,6)*pow(k,6)*pow(X1q,3)) ;
}

double J44 (const double & k, const double & X1q, const double & f1) {

	return  pow(E,-0.16666666666666666*f1*(2 + f1)*X1q*pow(k,2))*
     (0.1125541125541143 - 0.04812964812964716*f1*X1q*pow(k,2) + 
       X1q*pow(f1,2)*pow(k,2)*(-0.02406482406482358 + 0.014000814000813875*X1q*pow(k,2)) + 
       pow(f1,4)*(0.003500203500203469 - 0.004195078704882893*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) + 
       pow(f1,3)*(0.014000814000813875 - 0.002796719136588688*X1q*pow(k,2))*pow(k,4)*pow(X1q,2) - 
       0.0020975393524414465*pow(f1,5)*pow(k,6)*pow(X1q,3) - 0.000349589892073586*pow(f1,6)*pow(k,6)*pow(X1q,3)) ;
}


