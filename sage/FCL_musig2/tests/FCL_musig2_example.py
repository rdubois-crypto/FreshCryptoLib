##*************************************************************************************/
##/* Copyright (C) 2022 - Renaud Dubois - This file is part of cairo_musig2 project	 */
##/* License: This software is licensed under a dual BSD and GPL v2 license. 	 */
##/* See LICENSE file at the root folder of the project.				 */
##/* FILE: test_musig2.sage						             	  */
##/* 											  */
##/* 											  */
##/* DESCRIPTION: 2 round_multisignature Musig2 signatures 				  */
##/* This is a high level simulation for validation purpose				  */
##/* This file emulates a n-users aggregation and signatures, n being an input of sage */

##/* https:##eprint.iacr.org/2020/1261.pdf             				  */
##**************************************************************************************/

#note that the program expects the following value to be defined (currently done by makefile)
#to execute independtly, uncomment the following line:
#_MU=4;nb_users=3; size_message=2;seed=5;

from time import time
import sage.all
import subprocess
from sage.misc.sage_ostools import redirection


from FCL_musig2.FCL_musig2 import *

from FCL_musig2.FCL_musig2 import _MU
from FCL_common.FCL_io_printC import *

C_filepath='test_vector_musig2.c';
filep = open(C_filepath,'a');
_NTESTS=2;

_VERBOSE=false;

def print_verbose(x):
 if _VERBOSE:
    print(x)
    
def print_verbose_comment(comment, x):
 if _VERBOSE:
    print(comment, x)
    

def DisplaySave(comment, varname, counter, var):
	print_verbose_comment(comment, var);
	name=varname+str(counter);
	print_verbose_comment(name, hex(var));
	fprint_c_MSB(filep, "\n /*"+name+"*/",name, int(var),8,"");

		

with open(C_filepath, 'a') as file_out:
  
	print("   Test Musig2 with nb_users=",nb_users);
	filep.write("/* File automatically generated using test_musig2_example.sage, the  : ");
	with redirection(sys.stdout, filep):
		subprocess.call('date');
	filep = open(C_filepath,'a');

	print_verbose("\n\n********************************************************* \n*******************SAGEMATH:Simulation of a full Sign/Verif of Musig2:\n");
	print_verbose_comment("Simulation for a set of users of size:", nb_users);

	##***********************IO functions**********************************************/
	print_verbose_comment("\n*******************Generating Message of size",size_message);
	Fq=GF(Stark_order);
	message= [Fq(0) for i in range(5)] ;
	for i in range(0,size_message):
		message[i]=int(Fq.random_element());
	print_verbose(message);


	##***********************Key Generation and Aggregate******************************/

	print_verbose("\n*******************Generating Keys:\n");
	L=[];
	secrets=[];
	for i in range(0,nb_users):
		print_verbose_comment("\n***** user",i);
		[x,P]=Musig2_KeyGen(Curve, curve_Generator, Stark_order);#Concatenation of public keys
		
		DisplaySave("/*Secret key user", 'SecretKey_', i, x)
		
		DisplaySave("/*Public key ", 'PublicKey_X', i, int(P[0]))
		
		DisplaySave("/*", 'PublicKey_Y', i, int(P[1]))
		
		
#		print("/*Public key user",i,":\n",P,"*/");
#		name='PubX_'+str(i);
#		print(name, hex(var));
#		name='PubY_'+str(i);
#		print(name, hex(var));
		
		
		L=L+[int(P[0]),int(P[1])];	
		secrets=secrets+[x]; ##of course private keys are kept secret by users(simulation convenience)
		
	print_verbose_comment("\n /***Serialized public keys L:\n",L);
	print_verbose_comment("\n /***Serialized secret keys :\n",secrets);


	print_verbose_comment("\n ***Computing Aggregating coeffs ai:\n",L);	
	#vec_a=range(0,nb_users);	
	vec_a= [0 for i in range(nb_users)] ;
	for i in range(0,nb_users):	
		vec_a[i]=H_agg(L, nb_users, L[2*i], L[2*i+1], Stark_order);	
		#print("a[",i,"]=\n",vec_a[i]);
		print_verbose_comment("a[",i);
		print_verbose(vec_a[i]);
		
		name="a_"+str(i);
		fprint_c_MSB(filep, "\n /*"+name+"*/",name, int(vec_a[i]),8,"");
	

		
	print_verbose("\n*******************Aggregating Public Keys:\n");
	KeyAgg=Musig2_KeyAgg(Curve, curve_Generator, L, nb_users, Stark_order);
	print_verbose_comment("Aggregated Key:", KeyAgg);
	if(int(KeyAgg[1])&1==1): 
		print_verbose_comment("Wrong Key Agg, change Seed from", seed);
		exit();
	
	print_verbose("/*Aggregated Key: */");
	name='Key_Agg_X';
	#Print_C_MSB("\n /*"+name+"*/",name, int(KeyAgg[0]),8,"");
	fprint_c_MSB(filep, "\n /*"+name+"*/",name, int(KeyAgg[0]),8,"");
	
	name='Key_Agg_Y';
	#Print_C_MSB("\n /*"+name+"*/",name, int(KeyAgg[1]),8,"");
	fprint_c_MSB(filep, "\n /*"+name+"*/",name, int(KeyAgg[1]),8,"");
	

	##***********************Round 1 functions******************************/
	print_verbose("\n*******************Signature Round1:\n");
	infty_point=0*curve_Generator;

	vec_Rj=[0 for i in range(0,nb_users)];
	vec_nonces=[0 for i in range(0,nb_users)];	##of course nonces are kept secret by users(simulation convenience)


			
	for i in range(0,nb_users):#compute each users contribution
		[vec_nonces[i], vec_Rj[i]]=Musig2_Sign_Round1(Stark_order, nb_users,curve_Generator, _MU);	
		
	#Aggregate Round1 contributions	
	vec_R=Musig2_Sig1Agg(vec_Rj, nb_users, _MU);
		
	print_verbose_comment("All nonces:",vec_nonces);

	print_verbose_comment("Aggregated Signature=", vec_Rj);	

		
	##***********************Round 2 functions******************************/
	print_verbose("\n*******************Signature Round2:\n");

	vec_s=[0 for i in range(0,nb_users)];

	for i in range(0,nb_users):
		print_verbose_comment("\n***** User",i);
		[R,vec_s[i], c]=Musig2_Sign_Round2_all( 
			Curve, curve_Generator,Stark_order, nb_users, KeyAgg,#common parameters
			vec_a[i], secrets[i],					#user's data
			vec_R, vec_nonces[i],				#First round output
			message,  size_message, _MU);
			
		print_verbose_comment("s_: ", i);
		print_verbose(vec_s[i]);	
		#print("s_",i,":", vec_s[i]);
			
	s=Musig2_Sig2Agg(vec_s, Stark_order, nb_users);
			
	#print("Final Signature: \n R=", R,"\n s=", s,"\n c=", c);
	
	print_verbose_comment("Final Signature: \n R=", R);
	print_verbose_comment("s: \n R=", s);
	print_verbose_comment("c: \n R=", c);
	
	
	print_verbose_comment("With Key Agg:",KeyAgg);

	##***********************Verification functions******************************/
	print_verbose("\n*******************Verification :\n");
	print_verbose("/*Public KeyAgg */");
	print_verbose(hex(int(KeyAgg[0])));

	
	print_verbose("/*R part */");
	name='R_x';
	print_verbose(hex(int(R[0])));
	fprint_c_MSB(filep, "\n /*"+name+"*/",name, int(R[0]),8,"");
	
	name='R_y';
	
	fprint_c_MSB(filep, "\n /*"+name+"*/",name, int(R[1]),8,"");
	
	print_verbose("/*s part */");
	print_verbose(hex(int(s)));
	
	name='s';
	fprint_c_MSB(filep, "\n /*"+name+"*/",name, int(s),8,"");
	
	print_verbose("/*c part */");
	print_verbose(hex(int(c)));
	name='c';
	fprint_c_MSB(filep, "\n /*"+name+"*/",name, int(c),8,"");
	

	res_verif=Musig_Verif_Core(Curve, curve_Generator,R,s,KeyAgg, int(c));
	if(res_verif):
	  print("   OK");
	else:
	  print("KO");
		
	filep.close();




