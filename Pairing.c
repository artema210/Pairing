#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>

#define TRUE 1
#define FALSE 0
#define WSize 4

mpz_t z;
mpz_t prime;
mpz_t trace;
mpz_t order;
mpz_t r;
mpz_t r4;
mpz_t r8;
mpz_t r24;
unsigned int ECCPara_b;//b : y^2=x^3+ax+b

struct Fp4{
	//modular polynomial x^4+x^3+x^2+x+1

	mpz_t tau[4];
};

struct Fp8{
	//modular polynomial x^2-(tau+3)
	struct Fp4 theta[2];
};

struct Fp24{
	//modular polynomial x^6+x^5+x^4+x^3+x^2+x+1
	struct Fp8 omega[3];
};

	//y^2=x^3+b
struct EFp{
	mpz_t x,y;
	int Inf;
};


struct EFp4{
	struct Fp4 x,y;
	int Inf;
};

struct EFp8{
	struct Fp8 x,y;
	int Inf;
};

struct EFp24{
	struct Fp24 x,y;
	int Inf;
};

//----------Fp functions--------------------------------------------
void Fp_Pow(mpz_t Ans, mpz_t A, mpz_t B);	//Ans <- A^B
void Fp_Inv(mpz_t Ans, mpz_t A);			//Ans <- A^(-1)
void Fp_Sqrt(mpz_t Ans, mpz_t A);			//ANs <- A^(1/2)
//----------Fp4 functions--------------------------------------------
void Fp4_Init(struct Fp4 *X); // vector initialization
void Fp4_Clear(struct Fp4 *X);
void Fp4_Set(struct Fp4 *X, mpz_t *str);
void Fp4_Set1(struct Fp4 *X);	//X <- 1
void Fp4_Copy(struct Fp4 *X, struct Fp4 *A);	//X <- A
void Fp4_Add(struct Fp4 *Ans, struct Fp4 *A, struct Fp4 *B);	//Ans <- A+B
void Fp4_Sub(struct Fp4 *Ans, struct Fp4 *A, struct Fp4 *B);	//Ans <- A-B
void Fp4_Mul(struct Fp4 *Ans, struct Fp4 *A, struct Fp4 *B);	//Ans <- A*B
void Fp4_Div(struct Fp4 *Ans, struct Fp4 *A, struct Fp4 *B);	//Ans <- A/B
void Fp4_Neg(struct Fp4 *Ans, struct Fp4 *X);	//Ans <- -X
void Fp4_Sqr(struct Fp4 *Ans, struct Fp4 *X);	//Ans <- X^2
void Fp4_Sqrt(struct Fp4 *Ans, struct Fp4 *X);	//Ans <- X^(1/2)
void Fp4_Inv(struct Fp4 *Ans, struct Fp4 *X);	//Ans <- X^(-1)
void Fp4_Pow(struct Fp4 *Ans, struct Fp4 *A, mpz_t B);	//Ans <- A^B
void Fp4_Frob(struct Fp4 *Ans, struct Fp4 *X);	//Ans <- Frobenius Mapping of X
void Fp4_Rand(struct Fp4 *A);
void Fp4_Show(struct Fp4 *A);
int Fp4_Cmp(struct Fp4 *A, struct Fp4 *B);
int Fp4_Cmp_mpz(struct Fp4 *A, mpz_t B);
int Fp4_Legendre(struct Fp4 *A);

//----------Fp8 functions--------------------------------------------
void Fp8_Init(struct Fp8 *X);	// vector initialization
void Fp8_Clear(struct Fp8 *X);
void Fp8_Set(struct Fp8 *X, struct Fp4 *A);
void Fp8_Set1(struct Fp8 *X);	//X <- 1
void Fp8_Copy(struct Fp8 *X, struct Fp8 *A);	//X <- A
void Fp8_Add(struct Fp8 *Ans, struct Fp8 *A, struct Fp8 *B);	//Ans <- A+B
void Fp8_Sub(struct Fp8 *Ans, struct Fp8 *A, struct Fp8 *B);	//Ans <- A-B
void Fp8_Mul(struct Fp8 *Ans, struct Fp8 *A, struct Fp8 *B);	//Ans <- A*B
void Fp8_Div(struct Fp8 *Ans, struct Fp8 *A, struct Fp8 *B);	//Ans <- A/B
void Fp8_Neg(struct Fp8 *Ans, struct Fp8 *X);	//Ans <- -X
void Fp8_Sqr(struct Fp8 *Ans, struct Fp8 *X);	//Ans <- X^2
void Fp8_Sqrt(struct Fp8 *Ans, struct Fp8 *X);	//Ans <- X^(1/2)
void Fp8_Inv(struct Fp8 *Ans, struct Fp8 *X); //Ans <- X^(-1)
void Fp8_Pow(struct Fp8 *Ans, struct Fp8 *A, mpz_t B);	//Ans <- A^B
void Fp8_Frob(struct Fp8 *Ans, struct Fp8 *X);//Ans <- Frobenius Mapping of X
void Fp8_Rand(struct Fp8 *A);
void Fp8_Show(struct Fp8 *A);
int Fp8_Cmp(struct Fp8 *A, struct Fp8 *B);
int Fp8_Cmp_mpz(struct Fp8 *A, mpz_t B);
int Fp8_Legendre(struct Fp8 *A);

//----------Fp24 functions--------------------------------------------
void Fp24_Init(struct Fp24 *X);	// vector initialization
void Fp24_Clear(struct Fp24 *X);
void Fp24_Set(struct Fp24 *X, struct Fp8 *A);
void Fp24_Set1(struct Fp24 *X);	//X <- 1
void Fp24_Copy(struct Fp24 *X, struct Fp24 *A);	//X <- A
void Fp24_Add(struct Fp24 *Ans, struct Fp24 *A, struct Fp24 *B);	//Ans <- A+B
void Fp24_Sub(struct Fp24 *Ans, struct Fp24 *A, struct Fp24 *B);	//Ans <- A-B
void Fp24_Mul(struct Fp24 *Ans, struct Fp24 *A, struct Fp24 *B);	//Ans <- A*B
void Fp24_Div(struct Fp24 *Ans, struct Fp24 *A, struct Fp24 *B);	//Ans <- A/B
void Fp24_Neg(struct Fp24 *Ans, struct Fp24 *X);	//Ans <- -X
void Fp24_Sqr(struct Fp24 *Ans, struct Fp24 *X);	//Ans <- X^2
void Fp24_Sqrt(struct Fp24 *Ans, struct Fp24 *X);	//Ans <- X^(1/2)
void Fp24_Inv(struct Fp24 *Ans, struct Fp24 *X);	//Ans <- X^(-1)
void Fp24_Pow(struct Fp24 *Ans, struct Fp24 *A, mpz_t B);	//Ans <- A^B
void Fp24_Frob(struct Fp24 *Ans, struct Fp24 *X);	//Anx <- Frobenius Mapping of X
void Fp24_Rand(struct Fp24 *A);
void Fp24_Show(struct Fp24 *A);
int Fp24_Cmp(struct Fp24 *A, struct Fp24 *B);
int Fp24_Cmp_mpz(struct Fp24 *A, mpz_t B);
int Fp24_Legendre(struct Fp24 *A);

void Fp24_Check();	//simple check. Calculate A's pow(order);

//----------EFp functions--------------------------------------------
void EFp_Init(struct EFp *A);	// vector initialization
void EFp_Copy(struct EFp *A, struct EFp *B);	//A <- B
void EFp_SetInf(struct EFp *A);	//A <- Infinity
void EFp_Rand(struct EFp *A);
void EFp_Show(struct EFp *A);
void EFp_ECA(struct EFp *Ans, struct EFp *P1, struct EFp *P2);//ANS=P1+P2
void EFp_ECD(struct EFp *Ans, struct EFp *P);//ANS=2*P
void EFp_SCM(struct EFp *Ans, struct EFp *P,mpz_t j);//ANS=[j]P
//int EFp_cmp(struct EFp *A,struct EFp *B);//A==B->0
//void EFp_random_set(struct EFp *ANS);//generate random rational point

//----------EFp4 functions--------------------------------------------
void EFp4_Init(struct EFp4 *A);	// vector initialization
void EFp4_Copy(struct EFp4 *A, struct EFp4 *B);	//A <- B
void EFp4_SetInf(struct EFp4 *A);	//A <- Infinity
void EFp4_Rand(struct EFp4 *A);
void EFp4_Show(struct EFp4 *A);
void EFp4_ECA(struct EFp4 *Ans, struct EFp4 *P1, struct EFp4 *P2);//ANS=P1+P2
void EFp4_ECD(struct EFp4 *Ans, struct EFp4 *P);//ANS=2*P
void EFp4_SCM(struct EFp4 *Ans, struct EFp4 *P,mpz_t j);//ANS=[j]P

//----------EFp8 functions--------------------------------------------
void EFp8_Init(struct EFp8 *A);	// vector initialization
void EFp8_Copy(struct EFp8 *A, struct EFp8 *B);	//A <- B
void EFp8_SetInf(struct EFp8 *A);	// A <- Infinity
void EFp8_Rand(struct EFp8 *A);
void EFp8_Show(struct EFp8 *A);
void EFp8_ECA(struct EFp8 *Ans, struct EFp8 *P1, struct EFp8 *P2);//ANS=P1+P2
void EFp8_ECD(struct EFp8 *Ans, struct EFp8 *P);//ANS=2*P
void EFp8_SCM(struct EFp8 *Ans, struct EFp8 *P,mpz_t j);//ANS=[j]P

//----------EFp24 functions--------------------------------------------
void EFp24_Init(struct EFp24 *A);	// vector initialization
void EFp24_Copy(struct EFp24 *A, struct EFp24 *B);	//A <- B
void EFp24_SetInf(struct EFp24 *A);	//A <- Infinity
void EFp24_Rand(struct EFp24 *A);
void EFp24_Show(struct EFp24 *A);
void EFp24_ECA(struct EFp24 *Ans, struct EFp24 *P1, struct EFp24 *P2);//Ans=P1+P2
void EFp24_ECD(struct EFp24 *Ans, struct EFp24 *P);//Ans=2*P
void EFp24_SCM(struct EFp24 *Ans, struct EFp24 *P,mpz_t j);//Ans=[j]P
void EFp24_Inv(struct EFp24 *Ans, struct EFp24 *P);//Ans=[#EFp24-1]P=Inverse
void EFp24_Frob(struct EFp24 *Ans, struct EFp24 *P);	//Ans = Frobenius Map of P
void EFp24_G1_Rand(struct EFp24 *Ans);
void EFp24_G2_Rand(struct EFp24 *Ans);

//-------------Pairing functions--------------------------------
void Miller_Algo(struct Fp24 *Ans, struct EFp24 *P, struct EFp24 *Q, mpz_t roop);
void FDBL(struct Fp24 *Ans, struct EFp24 *R, struct EFp24 *T, struct EFp24 *Q);
void FADD(struct Fp24 *Ans, struct EFp24 *R, struct EFp24 *T, struct EFp24 *P, struct EFp24 *Q);
void Final_Exp(struct Fp24 *Ans, struct Fp24 *A);
void Tate_Pairing(struct Fp24 *Ans, struct EFp24 *G1, struct EFp24 *G2);

//-------------Fp functions--------------------------------
void Fp_Inv(mpz_t Ans, mpz_t A){
	mpz_t Ans_Copy,j;
	mpz_init(Ans_Copy);
	mpz_init(j);
	mpz_set(Ans_Copy, A);
	mpz_sub_ui(j,prime,2);
	int i,r;
	r=(int)mpz_sizeinbase(j,2);
	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			mpz_mul(Ans_Copy,Ans_Copy,Ans_Copy);
			mpz_mod(Ans_Copy,Ans_Copy,prime);
			mpz_mul(Ans_Copy,Ans_Copy,A);
		}else{
		mpz_mul(Ans_Copy,Ans_Copy,Ans_Copy);
		mpz_mod(Ans_Copy,Ans_Copy,prime);
		}
	}
	mpz_set(Ans,Ans_Copy);

	mpz_clear(Ans_Copy);
	mpz_clear(j);
}

void Fp_Pow(mpz_t Ans, mpz_t A, mpz_t B){
	if(mpz_cmp_ui(B,0)==0){
		mpz_set_ui(Ans,1);
		return;
	}
	mpz_t Ans_Copy;
	int i;
	int r;//bit

	mpz_init(Ans_Copy);
	r= (int)mpz_sizeinbase(B,2);
	mpz_set(Ans_Copy,A);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(B,i)==1){
			mpz_mul(Ans_Copy,Ans_Copy,Ans_Copy);	
			mpz_mul(Ans_Copy,Ans_Copy,A);			//(A*2)*A
			mpz_mod(Ans_Copy,Ans_Copy,prime);
		}else{
			mpz_mul(Ans_Copy,Ans_Copy,Ans_Copy);	//A*2
			mpz_mod(Ans_Copy,Ans_Copy,prime);
		}
	}
	mpz_mod(Ans_Copy,Ans_Copy,prime);
	mpz_set(Ans,Ans_Copy);

	mpz_clear(Ans_Copy);
}

void Fp_Pow2(mpz_t Ans, mpz_t A, mpz_t j){
	if(mpz_cmp_ui(j,0)==0){
		mpz_set_ui(Ans,1);
		return;
	}
	mpz_t Table[16],Ans_Copy;
	int i,r,k,Num;
	for(i=0;i<16;i++){
		mpz_init(Table[i]);
	}
	mpz_init(Ans_Copy);
	
	mpz_set_ui(Table[0],1);
	for(i=1;i<16;i++){
		mpz_mul(Table[i],Table[i-1],A);
		mpz_mod(Table[i],Table[i],prime);
	}
	
	//Size addjustment
	r=(int)mpz_sizeinbase(j,2);
	r+=(WSize-r%WSize)%WSize;
	r--;
	
	Num=0;
	for(k=1;k<WSize;k++){
		Num+=mpz_tstbit(j,r);
		Num=Num<<1;
		r--;
	}
	Num+=mpz_tstbit(j,r);
	r--;
	mpz_set_ui(Ans_Copy,1);
	mpz_mul(Ans_Copy,Ans_Copy,Table[Num]);
	for(i=r;i>=0;i-=WSize){
		for(k=0;k<WSize;k++){
			mpz_mul(Ans_Copy,Ans_Copy,Ans_Copy);
			mpz_mod(Ans_Copy,Ans_Copy,prime);
		}
		Num=0;
		for(k=1;k<WSize;k++){
			Num+=mpz_tstbit(j,r);
			Num=Num<<1;
			r--;
		}
		Num+=mpz_tstbit(j,r);
		r--;
		mpz_mul(Ans_Copy,Ans_Copy,Table[Num]);
		mpz_mod(Ans_Copy,Ans_Copy,prime);
	}
	mpz_set(Ans,Ans_Copy);
	
	for(i=0;i<16;i++){
		mpz_clear(Table[i]);
	}
	mpz_clear(Ans_Copy);
}

void Fp_Sqrt(mpz_t Ans, mpz_t A){
	if(mpz_legendre(A,prime)!=1){
		printf("No Answer\n");
		return;
	}

	mpz_t b,e,n,q,r,t,x,y,tmp1,tmp2,base;
	unsigned int m,Rand;
	gmp_randstate_t state;
	gmp_randinit_default(state);

	mpz_init(b);
	mpz_init(e);
	mpz_init(n);
	mpz_init(q);
	mpz_init(r);
	mpz_init(t);
	mpz_init(x);
	mpz_init(y);
	mpz_init(tmp1);
	mpz_init(tmp2);
	mpz_init(base);
	
	mpz_set_ui(base,2);
	mpz_set(n,A);
	while(mpz_legendre(n,prime)!=-1){
		Rand=rand();
		gmp_randseed_ui(state,Rand);
		mpz_urandomm(n, state, prime);
	}
	mpz_set_ui(n,21);
	mpz_set_ui(e,0);
	mpz_sub_ui(q,prime,1);
	
	while(mpz_odd_p(q)==0){
		mpz_add_ui(e,e,1);
		mpz_div_ui(q,q,2);
	}

	mpz_powm(y,n,q,prime);
	mpz_set(r,e);
	mpz_sub_ui(q,q,1);
	mpz_div_ui(q,q,2);
	mpz_powm(x,A,q,prime);
	mpz_set(tmp1,x);
	mpz_mul(x,x,A);
	mpz_mul(b,x,tmp1);
	mpz_mod(b,b,prime);

	while(mpz_cmp_ui(b,1)){
		m=-1;
		mpz_init(tmp2);
		while(mpz_cmp_ui(tmp2,1)){
			m++;
			mpz_ui_pow_ui(tmp1,2,m);
			mpz_powm(tmp2,b,tmp1,prime);
		}
		mpz_sub_ui(tmp1,r,m);
		mpz_sub_ui(tmp1,tmp1,1);
		mpz_powm(tmp2,base,tmp1,prime);
		mpz_powm(t,y,tmp2,prime);
		mpz_powm_ui(y,t,2,prime);
		mpz_set_ui(r,m);
		mpz_mul(x,x,t);
		mpz_mod(x,x,prime);
		mpz_mul(b,b,y);
		mpz_mod(b,b,prime);
	}
	
	mpz_set(Ans,x);
	mpz_mod(Ans,Ans,prime);
	
	mpz_clear(b);
	mpz_clear(e);
	mpz_clear(n);
	mpz_clear(q);
	mpz_clear(r);
	mpz_clear(t);
	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(tmp2);
	mpz_clear(tmp1);
	mpz_clear(base);

}

//-------------Fp4 functions--------------------------------
void Fp4_Init(struct Fp4 *X){
	int i;
	for(i=0;i<4;i++){ mpz_init(X->tau[i]); }
}

void Fp4_Clear(struct Fp4 *X){
	int i;
	for(i=0;i<4;i++){
		mpz_clear(X->tau[i]);
	}
}

void Fp4_Set(struct Fp4 *X, mpz_t *str){
	int i;
	for(i=0;i<4;i++){
		mpz_sub(X->tau[i], prime, str[i]);
	}
}

void Fp4_Set1(struct Fp4 *X){
	mpz_sub_ui(X->tau[0],prime,1);
	mpz_sub_ui(X->tau[1],prime,1);
	mpz_sub_ui(X->tau[2],prime,1);
	mpz_sub_ui(X->tau[3],prime,1);
}

void Fp4_Set_ui(struct Fp4 *A, unsigned int b){
	mpz_sub_ui(A->tau[0],prime,b);
	mpz_sub_ui(A->tau[1],prime,b);
	mpz_sub_ui(A->tau[2],prime,b);
	mpz_sub_ui(A->tau[3],prime,b);
}

/*
void Fp4_Set(struct Fp4 *A, struct Fp4 *B){
	int i;
	for(int i=0;i<4;i++){
		mpz_set_ui(X->tau[i],B->tau[i]);
	}
}
*/
void Fp4_Copy(struct Fp4 *X, struct Fp4 *A){
	int i;
	for(i=0;i<4;i++){
		mpz_set(X->tau[i],A->tau[i]);
	}
}

void Fp4_Add(struct Fp4 *Ans, struct Fp4 *A, struct Fp4 *B){
	int i;
	for(i=0;i<4;i++){
		mpz_add(Ans->tau[i],A->tau[i],B->tau[i]);
		mpz_fdiv_r(Ans->tau[i],Ans->tau[i],prime);
	}
}

void Fp4_Sub(struct Fp4 *Ans, struct Fp4 *A, struct Fp4 *B){
	int i;
	for(i=0;i<4;i++){
		mpz_sub(Ans->tau[i],A->tau[i],B->tau[i]);
		mpz_fdiv_r(Ans->tau[i],Ans->tau[i],prime);
	}
}

void Fp4_Mul(struct Fp4 *Ans, struct Fp4 *A, struct Fp4 *B){
	
	mpz_t tmp1,tmp2;
	struct Fp4 Ans_Copy;

	mpz_init(tmp1);
	mpz_init(tmp2);
	Fp4_Init(&Ans_Copy);

	mpz_sub(tmp1,A->tau[1],A->tau[3]);
	mpz_sub(tmp2,B->tau[1],B->tau[3]);
	mpz_mul(Ans_Copy.tau[0],tmp1,tmp2);
	mpz_sub(tmp1,A->tau[2],A->tau[3]);
	mpz_sub(tmp2,B->tau[2],B->tau[3]);
	mpz_mul(Ans_Copy.tau[1],tmp1,tmp2);
	mpz_sub(tmp1,A->tau[0],A->tau[1]);
	mpz_sub(tmp2,B->tau[0],B->tau[1]);
	mpz_mul(Ans_Copy.tau[2],tmp1,tmp2);
	mpz_sub(tmp1,A->tau[0],A->tau[2]);
	mpz_sub(tmp2,B->tau[0],B->tau[2]);
	mpz_mul(Ans_Copy.tau[3],tmp1,tmp2);
	
	mpz_add(tmp1,A->tau[0],A->tau[1]);
	mpz_sub(tmp1,tmp1,A->tau[2]);
	mpz_sub(tmp1,tmp1,A->tau[3]);
	mpz_add(tmp2,B->tau[0],B->tau[1]);
	mpz_sub(tmp2,tmp2,B->tau[2]);
	mpz_sub(tmp2,tmp2,B->tau[3]);
	mpz_mul(tmp1,tmp1,tmp2);

	mpz_sub(tmp1,tmp1,Ans_Copy.tau[0]);
	mpz_add(tmp1,tmp1,Ans_Copy.tau[1]);
	mpz_add(tmp1,tmp1,Ans_Copy.tau[2]);
	mpz_sub(tmp1,tmp1,Ans_Copy.tau[3]);

	mpz_sub(Ans_Copy.tau[0],tmp1,Ans_Copy.tau[0]);
	mpz_sub(Ans_Copy.tau[1],tmp1,Ans_Copy.tau[1]);
	mpz_sub(Ans_Copy.tau[2],tmp1,Ans_Copy.tau[2]);
	mpz_sub(Ans_Copy.tau[3],tmp1,Ans_Copy.tau[3]);

	mpz_mul(tmp1,A->tau[0],B->tau[0]);
	mpz_sub(Ans_Copy.tau[0],Ans_Copy.tau[0],tmp1);
	mpz_mul(tmp1,A->tau[1],B->tau[1]);
	mpz_sub(Ans_Copy.tau[1],Ans_Copy.tau[1],tmp1);
	mpz_mul(tmp1,A->tau[2],B->tau[2]);
	mpz_sub(Ans_Copy.tau[2],Ans_Copy.tau[2],tmp1);
	mpz_mul(tmp1,A->tau[3],B->tau[3]);
	mpz_sub(Ans_Copy.tau[3],Ans_Copy.tau[3],tmp1);
	
	mpz_fdiv_r(Ans_Copy.tau[0],Ans_Copy.tau[0],prime);
	mpz_fdiv_r(Ans_Copy.tau[1],Ans_Copy.tau[1],prime);
	mpz_fdiv_r(Ans_Copy.tau[2],Ans_Copy.tau[2],prime);
	mpz_fdiv_r(Ans_Copy.tau[3],Ans_Copy.tau[3],prime);

	Fp4_Copy(Ans,&Ans_Copy);

	mpz_clear(tmp1);
	mpz_clear(tmp2);
	Fp4_Clear(&Ans_Copy);
}

void Fp4_Div(struct Fp4 *Ans, struct Fp4 *A, struct Fp4 *B){
	struct Fp4 tmp;
	Fp4_Init(&tmp);

	Fp4_Inv(&tmp,B);
	Fp4_Mul(Ans,&tmp,A);
	
	Fp4_Clear(&tmp);
}

void Fp4_Neg(struct Fp4 *Ans,struct Fp4 *A){
	struct Fp4 Zero;
	Fp4_Init(&Zero);

	Fp4_Sub(Ans,&Zero,A);

	Fp4_Clear(&Zero);
}

void Fp4_Sqr(struct Fp4 *Ans, struct Fp4 *X){
	mpz_t tmp1,tmp2,tmp3;
	struct Fp4 Ans_Copy;

	mpz_init(tmp1);
	mpz_init(tmp2);
	mpz_init(tmp3);
	Fp4_Init(&Ans_Copy);

	mpz_sub(tmp1,X->tau[0],X->tau[1]);
	mpz_sub(tmp2,X->tau[2],X->tau[3]);
	mpz_mul(tmp1,tmp1,tmp2);
	mpz_mul_ui(tmp1,tmp1,2);
	mpz_sub(tmp2,X->tau[0],X->tau[2]);
	mpz_sub(tmp3,X->tau[1],X->tau[3]);
	mpz_mul(tmp2,tmp2,tmp3);
	mpz_mul_ui(tmp2,tmp2,2);
	
	mpz_mul_ui(tmp3,X->tau[0],2);
	mpz_sub(Ans_Copy.tau[0],tmp3,X->tau[2]);
	mpz_mul(Ans_Copy.tau[0],Ans_Copy.tau[0],X->tau[2]);
	mpz_mul_ui(tmp3,X->tau[1],2);
	mpz_sub(Ans_Copy.tau[1],tmp3,X->tau[0]);
	mpz_mul(Ans_Copy.tau[1],Ans_Copy.tau[1],X->tau[0]);
	mpz_mul_ui(tmp3,X->tau[2],2);
	mpz_sub(Ans_Copy.tau[2],tmp3,X->tau[3]);
	mpz_mul(Ans_Copy.tau[2],Ans_Copy.tau[2],X->tau[3]);
	mpz_mul_ui(tmp3,X->tau[3],2);
	mpz_sub(Ans_Copy.tau[3],tmp3,X->tau[1]);
	mpz_mul(Ans_Copy.tau[3],Ans_Copy.tau[3],X->tau[1]);

	mpz_sub(Ans_Copy.tau[0],tmp1,Ans_Copy.tau[0]);
	mpz_sub(Ans_Copy.tau[1],tmp2,Ans_Copy.tau[1]);
	mpz_sub(Ans_Copy.tau[2],tmp2,Ans_Copy.tau[2]);
	mpz_sub(Ans_Copy.tau[3],tmp1,Ans_Copy.tau[3]);

	mpz_fdiv_r(Ans_Copy.tau[0],Ans_Copy.tau[0],prime);
	mpz_fdiv_r(Ans_Copy.tau[1],Ans_Copy.tau[1],prime);
	mpz_fdiv_r(Ans_Copy.tau[2],Ans_Copy.tau[2],prime);
	mpz_fdiv_r(Ans_Copy.tau[3],Ans_Copy.tau[3],prime);

	Fp4_Copy(Ans,&Ans_Copy);

	mpz_clear(tmp1);
	mpz_clear(tmp2);
	mpz_clear(tmp3);
	Fp4_Clear(&Ans_Copy);
}

void Fp4_Sqrt(struct Fp4 *Ans,struct Fp4 *A){
	
	if(Fp4_Legendre(A)!=1){
		printf("No Answer\n");
		return;
	}
	
	struct Fp4 b,n,t,x,y,tmp1;
	mpz_t e,q,r,tmp2,base,set1;
	unsigned int m;
	
	Fp4_Init(&b);
	Fp4_Init(&n);
	Fp4_Init(&t);
	Fp4_Init(&x);
	Fp4_Init(&y);
	Fp4_Init(&tmp1);
	mpz_init(e);
	mpz_init(q);
	mpz_init(r);
	mpz_init(tmp2);
	mpz_init(base);
	mpz_init(set1);

	mpz_set_ui(set1,1);
	mpz_set_ui(base,2);
	Fp4_Copy(&n,A);
	
	while(Fp4_Legendre(&n)!=-1){
		Fp4_Rand(&n);
	}

	mpz_set_ui(e,0);
	mpz_pow_ui(q,prime,4);
	mpz_sub_ui(q,q,1);
	
	while(mpz_odd_p(q)==0){
		mpz_add_ui(e,e,1);
		mpz_div_ui(q,q,2);
	}

	Fp4_Pow(&y,&n,q);
	mpz_set(r,e);
	mpz_sub_ui(q,q,1);
	mpz_div_ui(q,q,2);
	Fp4_Pow(&x,A,q);
	Fp4_Copy(&tmp1,&x);
	Fp4_Mul(&x,&x,A);
	Fp4_Mul(&b,&x,&tmp1);
	Fp4_Init(&tmp1);

	while(Fp4_Cmp_mpz(&b,set1)){
		m=-1;
		Fp4_Init(&tmp1);
		while(Fp4_Cmp_mpz(&tmp1,set1)){
			m++;
			mpz_ui_pow_ui(tmp2,2,m);
			Fp4_Pow(&tmp1,&b,tmp2);
		}
		mpz_sub_ui(tmp2,r,m);
		mpz_sub_ui(tmp2,tmp2,1);
		Fp_Pow(tmp2,base,tmp2);
		Fp4_Pow(&t,&y,tmp2);
		Fp4_Mul(&y,&t,&t);
		mpz_set_ui(r,m);
		Fp4_Mul(&x,&x,&t);
		Fp4_Mul(&b,&b,&y);
	}
	
	Fp4_Copy(Ans,&x);
	
	Fp4_Clear(&b);
	Fp4_Clear(&n);
	Fp4_Clear(&t);
	Fp4_Clear(&x);
	Fp4_Clear(&y);
	Fp4_Clear(&tmp1);
	mpz_clear(e);
	mpz_clear(q);
	mpz_clear(r);
	mpz_clear(tmp2);
	mpz_clear(base);
	mpz_clear(set1);
}

void Fp4_Inv(struct Fp4 *Ans, struct Fp4 *X){
	mpz_t tmp;
	struct Fp4 tmp1,Ans_Copy;
	mpz_init(tmp);
	Fp4_Init(&tmp1);
	Fp4_Init(&Ans_Copy);
	
	Fp4_Frob(&tmp1,X);
	Fp4_Copy(&Ans_Copy,&tmp1);
	Fp4_Frob(&tmp1,&tmp1);
	Fp4_Mul(&Ans_Copy,&Ans_Copy,&tmp1);
	Fp4_Frob(&tmp1,&tmp1);
	Fp4_Mul(&Ans_Copy,&Ans_Copy,&tmp1);
	Fp4_Mul(&tmp1,&Ans_Copy,X);
	Fp_Inv(tmp,tmp1.tau[0]);
	
	mpz_set(tmp1.tau[0],tmp);
	mpz_set(tmp1.tau[1],tmp);
	mpz_set(tmp1.tau[2],tmp);
	mpz_set(tmp1.tau[3],tmp);

	Fp4_Mul(&Ans_Copy,&Ans_Copy,&tmp1);

	Fp4_Copy(Ans, &Ans_Copy);
	
	mpz_clear(tmp);
	Fp4_Clear(&tmp1);
	Fp4_Clear(&Ans_Copy);
}

void Fp4_Pow(struct Fp4 *Ans,struct Fp4 *A,mpz_t B){
	
	if(mpz_cmp_ui(B,0)==0){
		Fp4_Set1(Ans);
		return;
	}

	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(B,2);

	struct Fp4 Ans_Copy,tmp;
	Fp4_Init(&Ans_Copy);
	Fp4_Init(&tmp);
	Fp4_Copy(&Ans_Copy,A);
	Fp4_Copy(&tmp,A);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(B,i)==1){
			Fp4_Mul(&Ans_Copy,&Ans_Copy,&Ans_Copy);//a*2
			Fp4_Mul(&Ans_Copy,&Ans_Copy,&tmp);//*a
		}else{
			Fp4_Mul(&Ans_Copy,&Ans_Copy,&Ans_Copy);//a*2
		}
	}

	Fp4_Copy(Ans,&Ans_Copy);

	Fp4_Clear(&Ans_Copy);
	Fp4_Clear(&tmp);
}

void Fp4_Pow2(struct Fp4 *Ans, struct Fp4 *P, mpz_t j){
	struct Fp4 Table[16],Ans_Copy;
	int i,r,k,Num;
	for(i=0;i<16;i++){
		Fp4_Init(&Table[i]);
	}
	Fp4_Init(&Ans_Copy);
	
	Fp4_Set1(&Ans_Copy);
	Fp4_Set1(&Table[0]);
	for(i=1;i<16;i++){
		Fp4_Mul(&Table[i],&Table[i-1],P);
	}
	
	//Size addjustment
	r=(int)mpz_sizeinbase(j,2);
	r+=(WSize-r%WSize)%WSize;
	r--;
	
	Num=0;
	for(k=1;k<WSize;k++){
		Num+=mpz_tstbit(j,r);
		Num=Num<<1;
		r--;
	}
	Num+=mpz_tstbit(j,r);
	r--;
	Fp4_Mul(&Ans_Copy,&Ans_Copy,&Table[Num]);
	for(i=r;i>=0;i-=WSize){
		for(k=0;k<WSize;k++){
			Fp4_Sqr(&Ans_Copy,&Ans_Copy);
		}
		Num=0;
		for(k=1;k<WSize;k++){
			Num+=mpz_tstbit(j,r);
			Num=Num<<1;
			r--;
		}
		Num+=mpz_tstbit(j,r);
		r--;
		Fp4_Mul(&Ans_Copy,&Ans_Copy,&Table[Num]);
	}
	Fp4_Copy(Ans,&Ans_Copy);
	
	for(i=0;i<16;i++){
		Fp4_Clear(&Table[i]);
	}
	Fp4_Clear(&Ans_Copy);
}

void Fp4_Frob(struct Fp4 *Ans, struct Fp4 *X){
	struct Fp4 Ans_Copy;
	Fp4_Init(&Ans_Copy);
	mpz_set(Ans_Copy.tau[0],X->tau[2]);
	mpz_set(Ans_Copy.tau[1],X->tau[0]);
	mpz_set(Ans_Copy.tau[2],X->tau[3]);
	mpz_set(Ans_Copy.tau[3],X->tau[1]);
	Fp4_Copy(Ans,&Ans_Copy);

	Fp4_Clear(&Ans_Copy);
}

void Fp4_Rand(struct Fp4 *A){
	unsigned int i,r;
	gmp_randstate_t state;
	gmp_randinit_default(state);
	for(i=0;i<4;i++){
		r=rand();
		gmp_randseed_ui(state,r);
		mpz_urandomm(A->tau[i], state, prime);
	}
}

void Fp4_Show(struct Fp4 *A){
	int i;
	printf("(");
	for(i=0;i<4;i++){
		mpz_out_str(stdout, 10, A->tau[i]);
		printf(",");
	}
	printf(")\n");
}

int Fp4_Cmp(struct Fp4 *A,struct Fp4 *B){
	if(mpz_cmp(A->tau[0],B->tau[0])==0 && mpz_cmp(A->tau[1],B->tau[1])==0 && mpz_cmp(A->tau[2],B->tau[2])==0 && mpz_cmp(A->tau[3],B->tau[3])==0){
		return 0;
	}
	return 1;
}

int Fp4_Cmp_mpz(struct Fp4 *A,mpz_t B){
	mpz_t tmp;
	mpz_init(tmp);

	mpz_sub(tmp,prime,B);
	if(mpz_cmp(A->tau[0],tmp)==0 && mpz_cmp(A->tau[1],tmp)==0 && mpz_cmp(A->tau[2],tmp)==0 && mpz_cmp(A->tau[3],tmp)==0){
		mpz_clear(tmp);
		return 0;
	}
	mpz_clear(tmp);
	return 1;
}

int Fp4_Legendre(struct Fp4 *A){
	mpz_t i,cmp;
	struct Fp4 tmp;
	mpz_init(i);
	mpz_init(cmp);
	Fp4_Init(&tmp);

	mpz_set_ui(cmp,1);
	mpz_pow_ui(i,prime,4);
	mpz_sub_ui(i,i,1);
	mpz_tdiv_q_ui(i,i,2);
	Fp4_Pow(&tmp,A,i);

	if((Fp4_Cmp_mpz(&tmp,cmp))==0){
		Fp4_Clear(&tmp);
		mpz_clear(i);
		mpz_clear(cmp);
		return 1;
	}else{
		Fp4_Clear(&tmp);
		mpz_clear(i);
		mpz_clear(cmp);
		return -1;
	}
}

//-------------Fp8 functions--------------------------------

void Fp8_Init(struct Fp8 *X){
	Fp4_Init(&X->theta[0]);
	Fp4_Init(&X->theta[1]);
}

void Fp8_Clear(struct Fp8 *X){
	Fp4_Clear(&X->theta[0]);
	Fp4_Clear(&X->theta[1]);
}

void Fp8_Set(struct Fp8 *A, struct Fp4 *B){
	Fp4_Set(&A->theta[0], B[0].tau);
	Fp4_Set(&A->theta[1], B[1].tau);
}

void Fp8_Set1(struct Fp8 *X){
	Fp4_Set1(&X->theta[0]);
}

void Fp8_Set_ui(struct Fp8 *A,unsigned int b){
	Fp4_Set_ui(&A->theta[0],b);
	Fp4_Init(&A->theta[1]);
}

void Fp8_Copy(struct Fp8 *X, struct Fp8 *A){
	Fp4_Copy(&X->theta[0],&A->theta[0]);
	Fp4_Copy(&X->theta[1],&A->theta[1]);
}

void Fp8_Add(struct Fp8 *Ans, struct Fp8 *A, struct Fp8 *B){
	Fp4_Add(&Ans->theta[0], &A->theta[0], &B->theta[0]);
	Fp4_Add(&Ans->theta[1], &A->theta[1], &B->theta[1]);
}

void Fp8_Sub(struct Fp8 *Ans, struct Fp8 *A, struct Fp8 *B){
	Fp4_Sub(&Ans->theta[0], &A->theta[0], &B->theta[0]);
	Fp4_Sub(&Ans->theta[1], &A->theta[1], &B->theta[1]);
}

void Fp8_Mul(struct Fp8 *Ans, struct Fp8 *A, struct Fp8 *B){
	//Using Karatsuba Method
	struct Fp4 Theta,tmp1,tmp2;
	struct Fp8 Ans_Copy;
	Fp4_Init(&Theta);
	Fp4_Init(&tmp1);
	Fp4_Init(&tmp2);
	Fp8_Init(&Ans_Copy);
	mpz_sub_ui(Theta.tau[0], prime, 2);
	mpz_sub_ui(Theta.tau[1], prime, 3);
	mpz_sub_ui(Theta.tau[2], prime, 3);
	mpz_sub_ui(Theta.tau[3], prime, 3);
	/*
	mpz_set_str(Theta.tau[0], "25522886250087214",10);
	mpz_set_str(Theta.tau[1], "25522886250087213",10);
	mpz_set_str(Theta.tau[2], "25522886250087213",10);
	mpz_set_str(Theta.tau[3], "25522886250087213",10);
*/

	Fp4_Add(&tmp1, &A->theta[0], &A->theta[1]);
	Fp4_Add(&tmp2, &B->theta[0], &B->theta[1]);
	Fp4_Mul(&Ans_Copy.theta[1], &tmp1, &tmp2);

	Fp4_Mul(&tmp1, &A->theta[0], &B->theta[0]);
	Fp4_Mul(&tmp2, &A->theta[1], &B->theta[1]);
	Fp4_Sub(&Ans_Copy.theta[1], &Ans_Copy.theta[1], &tmp1);
	Fp4_Sub(&Ans_Copy.theta[1], &Ans_Copy.theta[1], &tmp2);

	Fp4_Mul(&Ans_Copy.theta[0], &Theta, &tmp2);
	Fp4_Add(&Ans_Copy.theta[0], &Ans_Copy.theta[0], &tmp1);

	Fp8_Copy(Ans,&Ans_Copy);

	Fp4_Clear(&Theta);
	Fp4_Clear(&tmp1);
	Fp4_Clear(&tmp2);
	Fp8_Clear(&Ans_Copy);
}

void Fp8_Div(struct Fp8 *Ans, struct Fp8 *A, struct Fp8 *B){
	struct Fp8 tmp;
	Fp8_Init(&tmp);

	Fp8_Inv(&tmp,B);
	Fp8_Mul(Ans,&tmp,A);
	
	Fp8_Clear(&tmp);
}

void Fp8_Neg(struct Fp8 *Ans,struct Fp8 *A){
	struct Fp8 Zero;
	Fp8_Init(&Zero);

	Fp8_Sub(Ans,&Zero,A);

	Fp8_Clear(&Zero);
}

void Fp8_Sqr(struct Fp8 *Ans, struct Fp8 *X){
	struct Fp4 Theta,tmp;
	struct Fp8 Ans_Copy;
	Fp4_Init(&Theta);
	Fp4_Init(&tmp);
	Fp8_Init(&Ans_Copy);
	mpz_sub_ui(Theta.tau[0], prime, 2);
	mpz_sub_ui(Theta.tau[1], prime, 3);
	mpz_sub_ui(Theta.tau[2], prime, 3);
	mpz_sub_ui(Theta.tau[3], prime, 3);
	
	Fp4_Add(&Ans_Copy.theta[1], &X->theta[0], &X->theta[1]);
	Fp4_Sqr(&Ans_Copy.theta[1], &Ans_Copy.theta[1]);
	Fp4_Sqr(&tmp, &X->theta[0]);
	Fp4_Sub(&Ans_Copy.theta[1], &Ans_Copy.theta[1], &tmp);
	Fp4_Copy(&Ans_Copy.theta[0], &tmp);
	Fp4_Sqr(&tmp, &X->theta[1]);
	Fp4_Sub(&Ans_Copy.theta[1], &Ans_Copy.theta[1], &tmp);
	Fp4_Mul(&tmp, &tmp, &Theta);
	Fp4_Add(&Ans_Copy.theta[0], &Ans_Copy.theta[0], &tmp);

	Fp8_Copy(Ans, &Ans_Copy);

	Fp4_Clear(&Theta);
	Fp4_Clear(&tmp);
	Fp8_Clear(&Ans_Copy);
}

void Fp8_Sqrt(struct Fp8 *Ans,struct Fp8 *A){
	
	if(Fp8_Legendre(A)!=1){
		printf("No Answer\n");
		return;
	}
	
	struct Fp8 b,n,t,x,y,tmp1;
	mpz_t e,q,r,tmp2,base,set1;
	unsigned int m;
	
	Fp8_Init(&b);
	Fp8_Init(&n);
	Fp8_Init(&t);
	Fp8_Init(&x);
	Fp8_Init(&y);
	Fp8_Init(&tmp1);
	mpz_init(e);
	mpz_init(q);
	mpz_init(r);
	mpz_init(tmp2);
	mpz_init(base);
	mpz_init(set1);

	mpz_set_ui(set1,1);
	mpz_set_ui(base,2);
	Fp8_Copy(&n,A);
	
	while(Fp8_Legendre(&n)!=-1){
		Fp8_Rand(&n);
	}
	mpz_set_ui(e,0);
	mpz_pow_ui(q,prime,8);
	mpz_sub_ui(q,q,1);
	
	while(mpz_odd_p(q)==0){
		mpz_add_ui(e,e,1);
		mpz_div_ui(q,q,2);
	}

	Fp8_Pow(&y,&n,q);
	mpz_set(r,e);
	mpz_sub_ui(q,q,1);
	mpz_div_ui(q,q,2);
	Fp8_Pow(&x,A,q);
	Fp8_Copy(&tmp1,&x);
	Fp8_Mul(&x,&x,A);
	Fp8_Mul(&b,&x,&tmp1);

	while(Fp8_Cmp_mpz(&b,set1)){
		m=-1;
		Fp8_Init(&tmp1);
		while(Fp8_Cmp_mpz(&tmp1,set1)){
			m++;
			mpz_ui_pow_ui(tmp2,2,m);
			Fp8_Pow(&tmp1,&b,tmp2);
		}
		mpz_sub_ui(tmp2,r,m);
		mpz_sub_ui(tmp2,tmp2,1);
		mpz_powm(tmp2,base,tmp2,prime);
		Fp8_Pow(&t,&y,tmp2);
		Fp8_Mul(&y,&t,&t);
		mpz_set_ui(r,m);
		Fp8_Mul(&x,&x,&t);
		Fp8_Mul(&b,&b,&y);
	}
	
	Fp8_Copy(Ans,&x);
	
	Fp8_Clear(&b);
	Fp8_Clear(&n);
	Fp8_Clear(&t);
	Fp8_Clear(&x);
	Fp8_Clear(&y);
	Fp8_Clear(&tmp1);
	mpz_clear(e);
	mpz_clear(q);
	mpz_clear(r);
	mpz_clear(tmp2);
	mpz_clear(base);
	mpz_clear(set1);

}

void Fp8_Inv(struct Fp8 *Ans, struct Fp8 *X){
	struct Fp8 tmp1,Ans_Copy;
	Fp8_Init(&tmp1);
	Fp8_Init(&Ans_Copy);

	Fp8_Frob(&Ans_Copy,X);
	Fp8_Mul(&tmp1,& Ans_Copy,X);
	Fp4_Inv(&tmp1.theta[0],&tmp1.theta[0]);

	Fp4_Mul(&Ans_Copy.theta[0],&Ans_Copy.theta[0],&tmp1.theta[0]);
	Fp4_Mul(&Ans_Copy.theta[1],&Ans_Copy.theta[1],&tmp1.theta[0]);

	Fp8_Copy(Ans,&Ans_Copy);
	
	Fp8_Clear(&tmp1);
	Fp8_Clear(&Ans_Copy);
}

void Fp8_Pow(struct Fp8 *Ans,struct Fp8 *A,mpz_t B){
	
	if(mpz_cmp_ui(B,0)==0){
		Fp8_Set1(Ans);
		return;
	}

	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(B,2);

	struct Fp8 Ans_Copy,tmp;
	Fp8_Init(&Ans_Copy);
	Fp8_Init(&tmp);
	Fp8_Copy(&Ans_Copy,A);
	Fp8_Copy(&tmp,A);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(B,i)==1){
			Fp8_Mul(&Ans_Copy,&Ans_Copy,&Ans_Copy);//a*2
			Fp8_Mul(&Ans_Copy,&Ans_Copy,&tmp);//*a
		}else{
			Fp8_Mul(&Ans_Copy,&Ans_Copy,&Ans_Copy);//a*2
		}
	}

	Fp8_Copy(Ans,&Ans_Copy);

	Fp8_Clear(&Ans_Copy);
	Fp8_Clear(&tmp);
}

void Fp8_Pow2(struct Fp8 *Ans, struct Fp8 *P, mpz_t j){
	struct Fp8 Table[16],Ans_Copy;
	int i,r,k,Num;
	for(i=0;i<16;i++){
		Fp8_Init(&Table[i]);
	}
	Fp8_Init(&Ans_Copy);
	
	Fp8_Set1(&Ans_Copy);
	Fp8_Set1(&Table[0]);
	for(i=1;i<16;i++){
		Fp8_Mul(&Table[i],&Table[i-1],P);
	}
	
	//Size addjustment
	r=(int)mpz_sizeinbase(j,2);
	r+=(WSize-r%WSize)%WSize;
	r--;
	
	Num=0;
	for(k=1;k<WSize;k++){
		Num+=mpz_tstbit(j,r);
		Num=Num<<1;
		r--;
	}
	Num+=mpz_tstbit(j,r);
	r--;
	Fp8_Mul(&Ans_Copy,&Ans_Copy,&Table[Num]);
	for(i=r;i>=0;i-=WSize){
		for(k=0;k<WSize;k++){
			Fp8_Sqr(&Ans_Copy,&Ans_Copy);
		}
		Num=0;
		for(k=1;k<WSize;k++){
			Num+=mpz_tstbit(j,r);
			Num=Num<<1;
			r--;
		}
		Num+=mpz_tstbit(j,r);
		r--;
		Fp8_Mul(&Ans_Copy,&Ans_Copy,&Table[Num]);
	}
	Fp8_Copy(Ans,&Ans_Copy);
	
	for(i=0;i<16;i++){
		Fp8_Clear(&Table[i]);
	}
	Fp8_Clear(&Ans_Copy);
}

void Fp8_Frob(struct Fp8 *Ans, struct Fp8 *X){
	
	struct Fp8 tmp;
	Fp8_Init(&tmp);

	Fp4_Sub(&tmp.theta[1],&tmp.theta[1],&X->theta[1]);
	mpz_set(Ans->theta[0].tau[0],X->theta[0].tau[0]);
	mpz_set(Ans->theta[0].tau[1],X->theta[0].tau[1]);
	mpz_set(Ans->theta[0].tau[2],X->theta[0].tau[2]);
	mpz_set(Ans->theta[0].tau[3],X->theta[0].tau[3]);
	mpz_set(Ans->theta[1].tau[0],tmp.theta[1].tau[0]);
	mpz_set(Ans->theta[1].tau[1],tmp.theta[1].tau[1]);
	mpz_set(Ans->theta[1].tau[2],tmp.theta[1].tau[2]);
	mpz_set(Ans->theta[1].tau[3],tmp.theta[1].tau[3]);

	Fp8_Clear(&tmp);
}

void Fp8_Rand(struct Fp8 *A){
	unsigned int i,j,r;
	gmp_randstate_t state;
	gmp_randinit_default(state);
	for(i=0;i<4;i++){
		for(j=0;j<2;j++){
			r=rand();
			gmp_randseed_ui(state,r);
			mpz_urandomm(A->theta[j].tau[i], state, prime);
		}
	}
}

void Fp8_Show(struct Fp8 *A){
	int i;
	printf("(");
	for(i=0;i<4;i++){
		mpz_out_str(stdout, 10, A->theta[0].tau[i]);
		printf(",");
	}
	printf(") + (");
	for(i=0;i<4;i++){
		mpz_out_str(stdout, 10, A->theta[1].tau[i]);
		printf(",");
	}
	printf(")Theta\n");
}

int Fp8_Cmp(struct Fp8 *A,struct Fp8 *B){
	int i;
	for(i=0;i<4;i++){
		if(mpz_cmp(A->theta[0].tau[i],B->theta[0].tau[i])||mpz_cmp(A->theta[1].tau[i],B->theta[1].tau[i])){
			return 1;
		}
	}
	return 0;
}

int Fp8_Cmp_mpz(struct Fp8 *A,mpz_t B){
	int i;
	mpz_t tmp;
	mpz_init(tmp);

	mpz_sub(tmp,prime,B);
	for(i=0;i<4;i++){
		if(mpz_cmp(A->theta[0].tau[i],tmp)||mpz_cmp_ui(A->theta[1].tau[i],0)){
			mpz_clear(tmp);
			return 1;
		}
	}
	mpz_clear(tmp);
	return 0;
}

int Fp8_Legendre(struct Fp8 *A){
	mpz_t i,cmp;
	struct Fp8 tmp;
	mpz_init(i);
	mpz_init(cmp);
	Fp8_Init(&tmp);

	mpz_set_ui(cmp,1);
	mpz_pow_ui(i,prime,8);
	mpz_sub_ui(i,i,1);
	mpz_tdiv_q_ui(i,i,2);
	Fp8_Pow(&tmp,A,i);

	if((Fp8_Cmp_mpz(&tmp,cmp))==0){
		Fp8_Clear(&tmp);
		mpz_clear(cmp);
		mpz_clear(i);
		return 1;
	}else{
		Fp8_Clear(&tmp);
		mpz_clear(i);
		mpz_clear(cmp);
		return -1;
	}

	return 0;
}

//-------------Fp24 functions--------------------------------

void Fp24_Init(struct Fp24 *X){
	Fp8_Init(&X->omega[0]);
	Fp8_Init(&X->omega[1]);
	Fp8_Init(&X->omega[2]);
}

void Fp24_Clear(struct Fp24 *X){
	Fp8_Clear(&X->omega[0]);
	Fp8_Clear(&X->omega[1]);
	Fp8_Clear(&X->omega[2]);
}

void Fp24_Set(struct Fp24 *A, struct Fp8 *B){
	Fp8_Set(&A->omega[0], B[0].theta);
	Fp8_Set(&A->omega[1], B[1].theta);
	Fp8_Set(&A->omega[2], B[2].theta);
}

void Fp24_Set1(struct Fp24 *X){
	mpz_set_ui(X->omega[0].theta[0].tau[0],1);
	mpz_set_ui(X->omega[0].theta[0].tau[1],1);
	mpz_set_ui(X->omega[0].theta[0].tau[2],1);
	mpz_set_ui(X->omega[0].theta[0].tau[3],1);
	mpz_set_ui(X->omega[1].theta[0].tau[0],1);
	mpz_set_ui(X->omega[1].theta[0].tau[1],1);
	mpz_set_ui(X->omega[1].theta[0].tau[2],1);
	mpz_set_ui(X->omega[1].theta[0].tau[3],1);
	mpz_set_ui(X->omega[2].theta[0].tau[0],1);
	mpz_set_ui(X->omega[2].theta[0].tau[1],1);
	mpz_set_ui(X->omega[2].theta[0].tau[2],1);
	mpz_set_ui(X->omega[2].theta[0].tau[3],1);
}

void Fp24_Set_ui(struct Fp24 *A, unsigned int b){
	mpz_set_ui(A->omega[0].theta[0].tau[0],b);
	mpz_set_ui(A->omega[0].theta[0].tau[1],b);
	mpz_set_ui(A->omega[0].theta[0].tau[2],b);
	mpz_set_ui(A->omega[0].theta[0].tau[3],b);
	mpz_set_ui(A->omega[1].theta[0].tau[0],b);
	mpz_set_ui(A->omega[1].theta[0].tau[1],b);
	mpz_set_ui(A->omega[1].theta[0].tau[2],b);
	mpz_set_ui(A->omega[1].theta[0].tau[3],b);
	mpz_set_ui(A->omega[2].theta[0].tau[0],b);
	mpz_set_ui(A->omega[2].theta[0].tau[1],b);
	mpz_set_ui(A->omega[2].theta[0].tau[2],b);
	mpz_set_ui(A->omega[2].theta[0].tau[3],b);
}

void Fp24_Copy(struct Fp24 *X, struct Fp24 *A){
	Fp8_Copy(&X->omega[0],&A->omega[0]);
	Fp8_Copy(&X->omega[1],&A->omega[1]);
	Fp8_Copy(&X->omega[2],&A->omega[2]);
}

void Fp24_Add(struct Fp24 *Ans, struct Fp24 *A, struct Fp24 *B){
	Fp8_Add(&Ans->omega[0], &A->omega[0], &B->omega[0]);
	Fp8_Add(&Ans->omega[1], &A->omega[1], &B->omega[1]);
	Fp8_Add(&Ans->omega[2], &A->omega[2], &B->omega[2]);
}

void Fp24_Sub(struct Fp24 *Ans, struct Fp24 *A, struct Fp24 *B){
	Fp8_Sub(&Ans->omega[0], &A->omega[0], &B->omega[0]);
	Fp8_Sub(&Ans->omega[1], &A->omega[1], &B->omega[1]);
	Fp8_Sub(&Ans->omega[2], &A->omega[2], &B->omega[2]);
}


/*
	//Multiplication by conventional method

void Fp24_Mul(struct Fp24 *Ans, struct Fp24 *A, struct Fp24 *B){
	//using Karatsuba Method

	struct Fp8 tmp,Theta;
	struct Fp24 Ans_Copy;

	Fp8_Init(&tmp);
	Fp8_Init(&Theta);
	Fp24_Init(&Ans_Copy);
	
	mpz_set_ui(Theta.theta[1].tau[0],p-1);
	mpz_set_ui(Theta.theta[1].tau[1],p-1);
	mpz_set_ui(Theta.theta[1].tau[2],p-1);
	mpz_set_ui(Theta.theta[1].tau[3],p-1);	//Theta <- (0,0,0,0)+(p-1,p-1,p-1,p-1)theta = theta
	
	Fp8_Add(&Ans_Copy.omega[0],&A->omega[1],&A->omega[2]);
	Fp8_Add(&tmp,&B->omega[1],&B->omega[2]);
	Fp8_Mul(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&tmp);
	Fp8_Add(&Ans_Copy.omega[1],&A->omega[0],&A->omega[1]);
	Fp8_Add(&tmp,&B->omega[0],&B->omega[1]);
	Fp8_Mul(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp);
	Fp8_Add(&Ans_Copy.omega[2],&A->omega[0],&A->omega[2]);
	Fp8_Add(&tmp,&B->omega[0],&B->omega[2]);
	Fp8_Mul(&Ans_Copy.omega[2],&Ans_Copy.omega[2],&tmp);
	
	Fp8_Mul(&tmp,&A->omega[1],&B->omega[1]);
	Fp8_Sub(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&tmp);
	Fp8_Sub(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp);
	Fp8_Add(&Ans_Copy.omega[2],&Ans_Copy.omega[2],&tmp);

	Fp8_Mul(&tmp,&A->omega[2],&B->omega[2]);
	Fp8_Sub(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&tmp);
	Fp8_Sub(&Ans_Copy.omega[2],&Ans_Copy.omega[2],&tmp);
	Fp8_Mul(&tmp,&tmp,&Theta);
	Fp8_Add(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp);

	Fp8_Mul(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&Theta);
	Fp8_Mul(&tmp,&A->omega[0],&B->omega[0]);
	Fp8_Add(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&tmp);
	Fp8_Sub(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp);
	Fp8_Sub(&Ans_Copy.omega[2],&Ans_Copy.omega[2],&tmp);

	Fp24_Copy(Ans, &Ans_Copy);
}
*/

void Fp24_Mul(struct Fp24 *Ans, struct Fp24 *A, struct Fp24 *B){
		//Multiplication by CVMA
	
	struct Fp8 tmp;
	struct Fp24 Ans_Copy;

	Fp8_Init(&tmp);
	Fp24_Init(&Ans_Copy);
	
	Fp8_Sub(&Ans_Copy.omega[0],&A->omega[0],&A->omega[1]);
	Fp8_Sub(&tmp,&B->omega[1],&B->omega[0]);
	Fp8_Mul(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&tmp);//Ans_Copy.omega[0] <- S1
	Fp8_Sub(&Ans_Copy.omega[1],&A->omega[1],&A->omega[2]);
	Fp8_Sub(&tmp,&B->omega[2],&B->omega[1]);
	Fp8_Mul(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp);//Ans_Copy.omega[1] <- S2
	Fp8_Sub(&Ans_Copy.omega[2],&A->omega[0],&A->omega[2]);
	Fp8_Sub(&tmp,&B->omega[2],&B->omega[0]);
	Fp8_Mul(&tmp,&Ans_Copy.omega[2],&tmp);//tmp <- S3

	Fp8_Add(&Ans_Copy.omega[2],&tmp,&Ans_Copy.omega[0]);
	Fp8_Add(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&Ans_Copy.omega[1]);
	Fp8_Add(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp);

	Fp8_Mul(&tmp,&A->omega[0],&B->omega[0]);
	Fp8_Sub(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&tmp);
	Fp8_Mul(&tmp,&A->omega[1],&B->omega[1]);
	Fp8_Sub(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp);
	Fp8_Mul(&tmp,&A->omega[2],&B->omega[2]);
	Fp8_Sub(&Ans_Copy.omega[2],&Ans_Copy.omega[2],&tmp);
	
	Fp24_Copy(Ans, &Ans_Copy);

	Fp8_Clear(&tmp);
	Fp24_Clear(&Ans_Copy);
}

void Fp24_Div(struct Fp24 *Ans, struct Fp24 *A, struct Fp24 *B){
	struct Fp24 tmp;
	Fp24_Init(&tmp);

	Fp24_Inv(&tmp,B);
	Fp24_Mul(Ans,&tmp,A);
	
	Fp24_Clear(&tmp);
}

void Fp24_Neg(struct Fp24 *Ans,struct Fp24 *A){
	struct Fp24 Zero;
	Fp24_Init(&Zero);

	Fp24_Sub(Ans,&Zero,A);

	Fp24_Clear(&Zero);
}

/*
void Fp24_Sqr(struct Fp24 *Ans, struct Fp24 *A){
	//using Karatsuba Method

	struct Fp8 tmp,Theta;
	struct Fp24 Ans_Copy;

	Fp8_Init(&tmp);
	Fp8_Init(&Theta);
	Fp24_Init(&Ans_Copy);

	mpz_set_ui(Theta.theta[1].tau[0],p-1);
	mpz_set_ui(Theta.theta[1].tau[1],p-1);
	mpz_set_ui(Theta.theta[1].tau[2],p-1);
	mpz_set_ui(Theta.theta[1].tau[3],p-1);	//Theta <- (0,0,0,0)+(p-1,p-1,p-1,p-1)theta = theta
	
	Fp8_Add(&Ans_Copy.omega[0],&A->omega[1],&A->omega[2]);
	Fp8_Sqr(&Ans_Copy.omega[0],&Ans_Copy.omega[0]);
	Fp8_Add(&Ans_Copy.omega[1],&A->omega[0],&A->omega[1]);
	Fp8_Sqr(&Ans_Copy.omega[1],&Ans_Copy.omega[1]);
	Fp8_Add(&Ans_Copy.omega[2],&A->omega[0],&A->omega[2]);
	Fp8_Sqr(&Ans_Copy.omega[2],&Ans_Copy.omega[2]);

	Fp8_Sqr(&tmp,&A->omega[1]);
	Fp8_Sub(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&tmp);
	Fp8_Sub(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp);
	Fp8_Add(&Ans_Copy.omega[2],&Ans_Copy.omega[2],&tmp);

	Fp8_Sqr(&tmp,&A->omega[2]);
	Fp8_Sub(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&tmp);
	Fp8_Sub(&Ans_Copy.omega[2],&Ans_Copy.omega[2],&tmp);
	Fp8_Mul(&tmp,&tmp,&Theta);
	Fp8_Add(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp);
	
	Fp8_Mul(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&Theta);
	Fp8_Sqr(&tmp,&A->omega[0]);
	Fp8_Add(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&tmp);
	Fp8_Sub(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp);
	Fp8_Sub(&Ans_Copy.omega[2],&Ans_Copy.omega[2],&tmp);
	Fp24_Copy(Ans, &Ans_Copy);
}
*/

void Fp24_Sqr(struct Fp24 *Ans, struct Fp24 *A){
	//using CVMA

	struct Fp8 tmp1, tmp2;
	struct Fp24 Ans_Copy;

	Fp8_Init(&tmp1);
	Fp8_Init(&tmp2);
	Fp24_Init(&Ans_Copy);

	Fp8_Sub(&Ans_Copy.omega[1],&A->omega[2],&A->omega[0]);
	Fp8_Sub(&tmp1,&A->omega[1],&A->omega[2]);
	Fp8_Mul(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp1);//Ans_Copy.omega[1] <-T3

	Fp8_Sub(&tmp1,&A->omega[1],&A->omega[0]);//tmp1 <-T1

	Fp8_Add(&tmp2,&A->omega[0],&A->omega[2]);
	Fp8_Mul(&tmp2,&tmp2,&tmp1);
	Fp8_Add(&Ans_Copy.omega[2],&tmp2,&Ans_Copy.omega[1]);//Ans_Copy.omega[2] <-T5

	Fp8_Add(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&Ans_Copy.omega[1]);//Ans_Copy.omega[1] <- 2T3

	Fp8_Add(&tmp2,&tmp1,&A->omega[2]);
	Fp8_Sqr(&tmp2,&tmp2);//tmp2 <-T2
	
	Fp8_Sqr(&Ans_Copy.omega[0],&A->omega[1]);
	Fp8_Sqr(&tmp1,&tmp1);
	Fp8_Add(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&tmp1);//Ans_Copy.omega[0] <-T4

	Fp8_Sub(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&Ans_Copy.omega[0]);//Ans_Copy.omega[1] <- 2T3-T4

	Fp8_Sub(&Ans_Copy.omega[0],&Ans_Copy.omega[2],&Ans_Copy.omega[0]);//Ans_Copy.omega[0] <- T5-T4
	Fp8_Sub(&Ans_Copy.omega[2],&Ans_Copy.omega[2],&tmp2);//Ans_Copy.omega[2] <- T5-T2
	
	Fp24_Copy(Ans,&Ans_Copy);
	
	Fp8_Clear(&tmp1);
	Fp8_Clear(&tmp2);
	Fp24_Clear(&Ans_Copy);
}

void Fp24_Sqrt(struct Fp24 *Ans,struct Fp24 *A){
	
	if(Fp24_Legendre(A)!=1){
		printf("No Answer\n");
		return;
	}
	
	struct Fp24 b,n,t,x,y,tmp1;
	mpz_t e,q,r,tmp2,base,set1;
	unsigned int m;
	
	Fp24_Init(&b);
	Fp24_Init(&n);
	Fp24_Init(&t);
	Fp24_Init(&x);
	Fp24_Init(&y);
	Fp24_Init(&tmp1);
	mpz_init(e);
	mpz_init(q);
	mpz_init(r);
	mpz_init(tmp2);
	mpz_init(base);
	mpz_init(set1);

	mpz_set_ui(set1,1);
	mpz_set_ui(base,2);
	Fp24_Copy(&n,A);
	
	while(Fp24_Legendre(&n)!=-1){
		Fp24_Rand(&n);
	}

	mpz_set_ui(e,0);
	mpz_pow_ui(q,prime,24);
	mpz_sub_ui(q,q,1);
	
	while(mpz_odd_p(q)==0){
		mpz_add_ui(e,e,1);
		mpz_div_ui(q,q,2);
	}

	Fp24_Pow(&y,&n,q);
	mpz_set(r,e);
	mpz_sub_ui(q,q,1);
	mpz_div_ui(q,q,2);
	Fp24_Pow(&x,A,q);
	Fp24_Copy(&tmp1,&x);
	Fp24_Mul(&x,&x,A);
	Fp24_Mul(&b,&x,&tmp1);
	Fp24_Init(&tmp1);

	while(Fp24_Cmp_mpz(&b,set1)){
		m=-1;
		Fp24_Init(&tmp1);
		while(Fp24_Cmp_mpz(&tmp1,set1)){
			m++;
			mpz_ui_pow_ui(tmp2,2,m);
			Fp24_Pow(&tmp1,&b,tmp2);
		}
		mpz_sub_ui(tmp2,r,m);
		mpz_sub_ui(tmp2,tmp2,1);
		Fp_Pow(tmp2,base,tmp2);
		Fp24_Pow(&t,&y,tmp2);
		Fp24_Mul(&y,&t,&t);
		mpz_set_ui(r,m);
		Fp24_Mul(&x,&x,&t);
		Fp24_Mul(&b,&b,&y);
	}
	
	Fp24_Copy(Ans,&x);
	
	Fp24_Clear(&b);
	Fp24_Clear(&n);
	Fp24_Clear(&t);
	Fp24_Clear(&x);
	Fp24_Clear(&y);
	Fp24_Clear(&tmp1);
	mpz_clear(e);
	mpz_clear(q);
	mpz_clear(r);
	mpz_clear(tmp2);
	mpz_clear(base);
	mpz_clear(set1);
}

/*
void Fp24_Inv(struct Fp24 *Ans, struct Fp24 *X){
	
	struct Fp24 Ans_Copy;
	char *str = "11110010110111100101101111011010001010100111100110010101100110001111110000000100100000111111001101010111111000011110101101010101010100111100111001010010111011010001001000100011011100101011011001101001000101110001111101010010001100000010111101111110011010101001000000010000111101011111111000011011111101001010010001011000110100011111010100101000010000000011000101001001001111110011100111100000111110110010011100011011001010011001101100010101000001101001011101110011100110010001100010101001000101111100001001101000011111111111011001001010010110101111111110000111111010001011000011000110111111010010011101101000001010010010111100001010111010010100111000110100101111101100000100101111110101100111011111111011011001100001101101101110001001110011101000011111";
	Fp24_Init(&Ans_Copy);
	mpz_set_str(Ans_Copy.omega[0].theta[0].tau[0],"25522886250087216",10);
	mpz_set_str(Ans_Copy.omega[0].theta[0].tau[1],"25522886250087216",10);
	mpz_set_str(Ans_Copy.omega[0].theta[0].tau[2],"25522886250087216",10);
	mpz_set_str(Ans_Copy.omega[0].theta[0].tau[3],"25522886250087216",10);
	unsigned int i;
	for(i=0;i<751;i++){
		if(str[i]=='1'){
			Fp24_Mul(&Ans_Copy,&Ans_Copy,X);
		}
		Fp24_Sqr(&Ans_Copy,&Ans_Copy);
	}
	Fp24_Mul(&Ans_Copy,&Ans_Copy,X);

	Fp24_Copy(Ans,&Ans_Copy);

}
*/

void Fp24_Inv(struct Fp24 *Ans, struct Fp24 *X){
	struct Fp24 tmp1,Ans_Copy;
	Fp24_Init(&tmp1);
	Fp24_Init(&Ans_Copy);

	Fp24_Frob(&Ans_Copy,X);
	Fp24_Frob(&tmp1,&Ans_Copy);
	Fp24_Mul(&Ans_Copy,&tmp1,&Ans_Copy);
	Fp24_Mul(&tmp1,&Ans_Copy,X);
	Fp8_Inv(&tmp1.omega[0],&tmp1.omega[0]);

	mpz_neg(tmp1.omega[0].theta[0].tau[0],tmp1.omega[0].theta[0].tau[0]);
	mpz_neg(tmp1.omega[0].theta[0].tau[1],tmp1.omega[0].theta[0].tau[1]);
	mpz_neg(tmp1.omega[0].theta[0].tau[2],tmp1.omega[0].theta[0].tau[2]);
	mpz_neg(tmp1.omega[0].theta[0].tau[3],tmp1.omega[0].theta[0].tau[3]);
	mpz_neg(tmp1.omega[0].theta[1].tau[0],tmp1.omega[0].theta[1].tau[0]);
	mpz_neg(tmp1.omega[0].theta[1].tau[1],tmp1.omega[0].theta[1].tau[1]);
	mpz_neg(tmp1.omega[0].theta[1].tau[2],tmp1.omega[0].theta[1].tau[2]);
	mpz_neg(tmp1.omega[0].theta[1].tau[3],tmp1.omega[0].theta[1].tau[3]);

	Fp8_Mul(&Ans_Copy.omega[0],&Ans_Copy.omega[0],&tmp1.omega[0]);
	Fp8_Mul(&Ans_Copy.omega[1],&Ans_Copy.omega[1],&tmp1.omega[0]);
	Fp8_Mul(&Ans_Copy.omega[2],&Ans_Copy.omega[2],&tmp1.omega[0]);

	Fp24_Copy(Ans,&Ans_Copy);
	
	Fp24_Clear(&tmp1);
	Fp24_Clear(&Ans_Copy);
}

void Fp24_Pow(struct Fp24 *Ans,struct Fp24 *A,mpz_t B){
	
	if(mpz_cmp_ui(B,0)==0){
		Fp24_Set1(Ans);
		return;
	}

	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(B,2);

	struct Fp24 Ans_Copy,tmp;
	Fp24_Init(&Ans_Copy);
	Fp24_Init(&tmp);
	Fp24_Copy(&Ans_Copy,A);
	Fp24_Copy(&tmp,A);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(B,i)==1){
			Fp24_Mul(&Ans_Copy,&Ans_Copy,&Ans_Copy);//a*2
			Fp24_Mul(&Ans_Copy,&Ans_Copy,&tmp);//*a
		}else{
			Fp24_Mul(&Ans_Copy,&Ans_Copy,&Ans_Copy);//a*2
		}
	}

	Fp24_Copy(Ans,&Ans_Copy);

	Fp24_Clear(&Ans_Copy);
	Fp24_Clear(&tmp);
}

void Fp24_Pow2(struct Fp24 *Ans, struct Fp24 *P, mpz_t j){
	struct Fp24 Table[16],Ans_Copy;
	int i,r,k,Num;
	for(i=0;i<16;i++){
		Fp24_Init(&Table[i]);
	}
	Fp24_Init(&Ans_Copy);
	
	Fp24_Set1(&Ans_Copy);
	Fp24_Set1(&Table[0]);
	for(i=1;i<16;i++){
		Fp24_Mul(&Table[i],&Table[i-1],P);
	}
	
	//Size addjustment
	r=(int)mpz_sizeinbase(j,2);
	r+=(WSize-r%WSize)%WSize;
	r--;
	
	Num=0;
	for(k=1;k<WSize;k++){
		Num+=mpz_tstbit(j,r);
		Num=Num<<1;
		r--;
	}
	Num+=mpz_tstbit(j,r);
	r--;
	Fp24_Mul(&Ans_Copy,&Ans_Copy,&Table[Num]);
	for(i=r;i>=0;i-=WSize){
		for(k=0;k<WSize;k++){
			Fp24_Sqr(&Ans_Copy,&Ans_Copy);
		}
		Num=0;
		for(k=1;k<WSize;k++){
			Num+=mpz_tstbit(j,r);
			Num=Num<<1;
			r--;
		}
		Num+=mpz_tstbit(j,r);
		r--;
		Fp24_Mul(&Ans_Copy,&Ans_Copy,&Table[Num]);
	}
	Fp24_Copy(Ans,&Ans_Copy);
	
	for(i=0;i<16;i++){
		Fp24_Clear(&Table[i]);
	}
	Fp24_Clear(&Ans_Copy);
}

void Fp24_Frob(struct Fp24 *Ans, struct Fp24 *X){
	struct Fp24 Ans_Copy;
	Fp24_Init(&Ans_Copy);

	mpz_set(Ans_Copy.omega[1].theta[0].tau[0],X->omega[0].theta[0].tau[0]);
	mpz_set(Ans_Copy.omega[1].theta[0].tau[1],X->omega[0].theta[0].tau[1]);
	mpz_set(Ans_Copy.omega[1].theta[0].tau[2],X->omega[0].theta[0].tau[2]);
	mpz_set(Ans_Copy.omega[1].theta[0].tau[3],X->omega[0].theta[0].tau[3]);
	mpz_set(Ans_Copy.omega[1].theta[1].tau[0],X->omega[0].theta[1].tau[0]);
	mpz_set(Ans_Copy.omega[1].theta[1].tau[1],X->omega[0].theta[1].tau[1]);
	mpz_set(Ans_Copy.omega[1].theta[1].tau[2],X->omega[0].theta[1].tau[2]);
	mpz_set(Ans_Copy.omega[1].theta[1].tau[3],X->omega[0].theta[1].tau[3]);
	mpz_set(Ans_Copy.omega[2].theta[0].tau[0],X->omega[1].theta[0].tau[0]);
	mpz_set(Ans_Copy.omega[2].theta[0].tau[1],X->omega[1].theta[0].tau[1]);
	mpz_set(Ans_Copy.omega[2].theta[0].tau[2],X->omega[1].theta[0].tau[2]);
	mpz_set(Ans_Copy.omega[2].theta[0].tau[3],X->omega[1].theta[0].tau[3]);
	mpz_set(Ans_Copy.omega[2].theta[1].tau[0],X->omega[1].theta[1].tau[0]);
	mpz_set(Ans_Copy.omega[2].theta[1].tau[1],X->omega[1].theta[1].tau[1]);
	mpz_set(Ans_Copy.omega[2].theta[1].tau[2],X->omega[1].theta[1].tau[2]);
	mpz_set(Ans_Copy.omega[2].theta[1].tau[3],X->omega[1].theta[1].tau[3]);
	mpz_set(Ans_Copy.omega[0].theta[0].tau[0],X->omega[2].theta[0].tau[0]);
	mpz_set(Ans_Copy.omega[0].theta[0].tau[1],X->omega[2].theta[0].tau[1]);
	mpz_set(Ans_Copy.omega[0].theta[0].tau[2],X->omega[2].theta[0].tau[2]);
	mpz_set(Ans_Copy.omega[0].theta[0].tau[3],X->omega[2].theta[0].tau[3]);
	mpz_set(Ans_Copy.omega[0].theta[1].tau[0],X->omega[2].theta[1].tau[0]);
	mpz_set(Ans_Copy.omega[0].theta[1].tau[1],X->omega[2].theta[1].tau[1]);
	mpz_set(Ans_Copy.omega[0].theta[1].tau[2],X->omega[2].theta[1].tau[2]);
	mpz_set(Ans_Copy.omega[0].theta[1].tau[3],X->omega[2].theta[1].tau[3]);

	Fp24_Copy(Ans,&Ans_Copy);
	Fp24_Clear(&Ans_Copy);
}

void Fp24_Frobp(struct Fp24 *Ans, struct Fp24 *X){
	struct Fp24 tmp;
	Fp24_Init(&tmp);

	Fp24_Pow(&tmp,X,prime);
	Fp24_Copy(Ans,&tmp);
	Fp24_Clear(&tmp);
}

void Fp24_Rand(struct Fp24 *A){
	unsigned int i,j,k,r;
	gmp_randstate_t state;
	gmp_randinit_default(state);
	for(i=0;i<4;i++){
		for(j=0;j<2;j++){
			for(k=0;k<3;k++){
				r=rand();
				gmp_randseed_ui(state,r);
				mpz_urandomm(A->omega[k].theta[j].tau[i], state, prime);
			}
		}
	}
}

void Fp24_Show(struct Fp24 *A){
	int i;
	printf("((");
	for(i=0;i<4;i++){
		mpz_out_str(stdout, 10, A->omega[0].theta[0].tau[i]);
		printf(",");
	}
	printf("),(");
	for(i=0;i<4;i++){
		mpz_out_str(stdout, 10, A->omega[0].theta[1].tau[i]);
		printf(",");
	}
	printf("))\n((");
	for(i=0;i<4;i++){
		mpz_out_str(stdout, 10, A->omega[1].theta[0].tau[i]);
		printf(",");
	}
	printf("),(");
	for(i=0;i<4;i++){
		mpz_out_str(stdout, 10, A->omega[1].theta[1].tau[i]);
		printf(",");
	}
	printf("))\n((");
	for(i=0;i<4;i++){
		mpz_out_str(stdout, 10, A->omega[2].theta[0].tau[i]);
		printf(",");
	}
	printf("),(");
	for(i=0;i<4;i++){
		mpz_out_str(stdout, 10, A->omega[2].theta[1].tau[i]);
		printf(",");
	}
	printf("))\n");
}

int Fp24_Cmp(struct Fp24 *A,struct Fp24 *B){
	int i,j;
	for(i=0;i<4;i++){
		for(j=0;j<3;j++){
			if(mpz_cmp(A->omega[j].theta[0].tau[i],B->omega[j].theta[0].tau[i])||mpz_cmp(A->omega[j].theta[1].tau[i],B->omega[j].theta[1].tau[i])){
				return 1;
			}
		}
	}
	return 0;
}

int Fp24_Cmp_mpz(struct Fp24 *A,mpz_t B){
	int i,j;
	
	for(i=0;i<4;i++){
		for(j=0;j<3;j++){
			if(mpz_cmp(A->omega[j].theta[0].tau[i],B)||mpz_cmp_ui(A->omega[j].theta[1].tau[i],0)){
				return 1;
			}
		}
	}
	return 0;
}

int Fp24_Legendre(struct Fp24 *A){
	mpz_t i,cmp;
	struct Fp24 tmp;
	mpz_init(i);
	mpz_init(cmp);
	Fp24_Init(&tmp);

	mpz_set_ui(cmp,1);
	mpz_pow_ui(i,prime,24);
	mpz_sub_ui(i,i,1);
	mpz_tdiv_q_ui(i,i,2);
	Fp24_Pow(&tmp,A,i);

	if((Fp24_Cmp_mpz(&tmp,cmp))==0){
		Fp24_Clear(&tmp);
		mpz_clear(i);
		mpz_clear(cmp);
		return 1;
	}else{
		Fp24_Clear(&tmp);
		mpz_clear(i);
		mpz_clear(cmp);
		return -1;
	}
	return 0;
}

//-------------EFp functions--------------------------------
void EFp_Init(struct EFp *A){
	mpz_init(A->x);
	mpz_init(A->y);
	A->Inf=FALSE;
}

void EFp_Copy(struct EFp *A, struct EFp *B){
	mpz_set(A->x,B->x);
	mpz_set(A->y,B->y);
	A->Inf=B->Inf;
}

void EFp_SetInf(struct EFp *A){
	mpz_set_ui(A->x,0);
	mpz_set_ui(A->y,0);
	A->Inf=TRUE;
}

void EFp_Clear(struct EFp *A){
	mpz_clear(A->x);
	mpz_clear(A->y);
}

void EFp_Rand(struct EFp *A){
	unsigned int r;
	struct EFp Ans_Copy;
	mpz_t tmp;
	EFp_Init(&Ans_Copy);
	mpz_init(tmp);
	gmp_randstate_t state;
	gmp_randinit_default(state);

	while(mpz_legendre(tmp,prime)!=1){
		r=rand();
		gmp_randseed_ui(state,r);
		mpz_urandomm(Ans_Copy.x,state,prime);
		mpz_powm_ui(tmp,Ans_Copy.x,3,prime);
		mpz_add_ui(tmp,tmp,ECCPara_b);
	}
	Fp_Sqrt(Ans_Copy.y,tmp);
	EFp_Copy(A,&Ans_Copy);

	mpz_clear(tmp);
	EFp_Clear(&Ans_Copy);
}

void EFp_Show(struct EFp *A){
	gmp_printf("(%Zd,%Zd)\n",A->x,A->y);
}

void EFp_ECA(struct EFp *Ans, struct EFp *P, struct EFp *Q){
	if(P->Inf==TRUE){
		EFp_Copy(Ans, Q);
		return;
	}else if(Q->Inf==TRUE){
		EFp_Copy(Ans, P);
		return;
	}else if(mpz_cmp(P->x,Q->x)==0&&mpz_cmp(P->y,Q->y)){
		EFp_SetInf(Ans);
		return;
	}else if(mpz_cmp(P->x,Q->x)==0&&mpz_cmp(P->y,Q->y)==0){
		EFp_ECD(Ans,P);
		return;
	}
	mpz_t lambda,tmp;
	struct EFp Ans_Copy;

	mpz_init(lambda);
	mpz_init(tmp);
	EFp_Init(&Ans_Copy);

	mpz_sub(tmp,P->x,Q->x);
	Fp_Inv(tmp,tmp);
	mpz_sub(lambda,P->y,Q->y);
	mpz_mul(lambda,lambda,tmp);
	mpz_mod(lambda,lambda,prime);

	mpz_mul(Ans_Copy.x,lambda,lambda);
	mpz_sub(Ans_Copy.x,Ans_Copy.x,P->x);
	mpz_sub(Ans_Copy.x,Ans_Copy.x,Q->x);
	mpz_mod(Ans_Copy.x,Ans_Copy.x,prime);

	mpz_sub(tmp,P->x,Ans_Copy.x);
	mpz_mul(Ans_Copy.y,lambda,tmp);
	mpz_sub(Ans_Copy.y,Ans_Copy.y,P->y);
	mpz_mod(Ans_Copy.y,Ans_Copy.y,prime);

	EFp_Copy(Ans,&Ans_Copy);

	mpz_clear(lambda);
	mpz_clear(tmp);
	EFp_Clear(&Ans_Copy);
}

void EFp_ECD(struct EFp *Ans, struct EFp *P){
	if(P->Inf==TRUE){
		EFp_SetInf(Ans);
		return;
	}else if(mpz_cmp_ui(P->y,0)==0){
		EFp_SetInf(Ans);
		return;
	}
	mpz_t lambda,tmp;
	struct EFp Ans_Copy;

	mpz_init(lambda);
	mpz_init(tmp);
	EFp_Init(&Ans_Copy);

	mpz_mul_ui(tmp,P->y,2);
	Fp_Inv(tmp,tmp);
	mpz_mul(lambda,P->x,P->x);
	mpz_mul_ui(lambda,lambda,3);
	mpz_mul(lambda,lambda,tmp);
	mpz_mod(lambda,lambda,prime);

	mpz_mul(Ans_Copy.x,lambda,lambda);
	mpz_mul_ui(tmp,P->x,2);
	mpz_sub(Ans_Copy.x,Ans_Copy.x,tmp);
	mpz_mod(Ans_Copy.x,Ans_Copy.x,prime);

	mpz_sub(tmp,P->x,Ans_Copy.x);
	mpz_mul(Ans_Copy.y,lambda,tmp);
	mpz_sub(Ans_Copy.y,Ans_Copy.y,P->y);
	mpz_mod(Ans_Copy.y,Ans_Copy.y,prime);

	EFp_Copy(Ans,&Ans_Copy);

	mpz_clear(lambda);
	mpz_clear(tmp);
	EFp_Clear(&Ans_Copy);
}

//Binary Method
void EFp_SCM(struct EFp *Ans, struct EFp *P, mpz_t j){
	int i;
	int r;
	r = (int)mpz_sizeinbase(j,2);
	struct EFp Ans_Copy;
	EFp_Init(&Ans_Copy);
	EFp_Copy(&Ans_Copy,P);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp_ECD(&Ans_Copy,&Ans_Copy);
			EFp_ECA(&Ans_Copy,&Ans_Copy,P);
		}else{
			EFp_ECD(&Ans_Copy,&Ans_Copy);
		}
	}

	EFp_Copy(Ans,&Ans_Copy);
	EFp_Clear(&Ans_Copy);
}

//Window Method
void EFp_SCM2(struct EFp *Ans, struct EFp *P, mpz_t j){
	int TabSize = 2<<(WSize-1);
	struct EFp Table[TabSize],Ans_Copy;
	int i,r,k,Num;
	for(i=0;i<TabSize;i++){
		EFp_Init(&Table[i]);
	}
	EFp_Init(&Ans_Copy);
	
	EFp_SetInf(&Ans_Copy);
	EFp_SetInf(&Table[0]);
	for(i=1;i<TabSize;i++){
		EFp_ECA(&Table[i],&Table[i-1],P);
	}
	
	//Size addjustment
	r=(int)mpz_sizeinbase(j,2);
	r+=(WSize-r%WSize)%WSize;
	r--;
	
	Num=0;
	for(k=1;k<WSize;k++){
		Num+=mpz_tstbit(j,r);
		Num=Num<<1;
		r--;
	}
	Num+=mpz_tstbit(j,r);
	r--;
	EFp_ECA(&Ans_Copy,&Ans_Copy,&Table[Num]);
	for(i=r;i>=0;i-=WSize){
		for(k=0;k<WSize;k++){
			EFp_ECD(&Ans_Copy,&Ans_Copy);
		}
		Num=0;
		for(k=1;k<WSize;k++){
			Num+=mpz_tstbit(j,r);
			Num=Num<<1;
			r--;
		}
		Num+=mpz_tstbit(j,r);
		r--;
		EFp_ECA(&Ans_Copy,&Ans_Copy,&Table[Num]);
	}
	EFp_Copy(Ans,&Ans_Copy);
	
	for(i=0;i<TabSize;i++){
		EFp_Clear(&Table[i]);
	}
	EFp_Clear(&Ans_Copy);
}

//-------------EFp4 functions--------------------------------
void EFp4_Init(struct EFp4 *A){
	Fp4_Init(&A->x);
	Fp4_Init(&A->y);
	A->Inf=FALSE;
}

void EFp4_Copy(struct EFp4 *A,struct EFp4 *B){
	Fp4_Copy(&A->x,&B->x);
	Fp4_Copy(&A->y,&B->y);
	A->Inf=B->Inf;
}

void EFp4_SetInf(struct EFp4 *A){
	Fp4_Init(&A->x);
	Fp4_Init(&A->y);
	A->Inf=TRUE;
}

void EFp4_Clear(struct EFp4 *A){
	Fp4_Clear(&A->x);
	Fp4_Clear(&A->y);
}

void EFp4_Rand(struct EFp4 *A){
	struct EFp4 Ans_Copy;
	struct Fp4 tmp,ECC_b;
	mpz_t j;
	EFp4_Init(&Ans_Copy);
	Fp4_Init(&tmp);
	Fp4_Init(&ECC_b);
	mpz_init(j);

	Fp4_Set_ui(&ECC_b,ECCPara_b);
	mpz_set_ui(j,3);
	do{
		Fp4_Rand(&Ans_Copy.x);
		Fp4_Pow(&tmp,&Ans_Copy.x,j);
		Fp4_Add(&tmp,&tmp,&ECC_b);
	}while(Fp4_Legendre(&tmp)!=1);

	Fp4_Sqrt(&Ans_Copy.y,&tmp);
	EFp4_Copy(A,&Ans_Copy);

	EFp4_Clear(&Ans_Copy);
	Fp4_Clear(&tmp);
	Fp4_Clear(&ECC_b);
	mpz_clear(j);
}

void EFp4_Show(struct EFp4 *A){
	gmp_printf("((%Zd,%Zd,%Zd,%Zd),(%Zd,%Zd,%Zd,%Zd))\n",A->x.tau[0],A->x.tau[1],A->x.tau[2],A->x.tau[3],A->y.tau[0],A->y.tau[1],A->y.tau[2],A->y.tau[3]);
}

void EFp4_ECA(struct EFp4 *Ans, struct EFp4 *P, struct EFp4 *Q){
	if(P->Inf==TRUE){
		EFp4_Copy(Ans, Q);
		return;
	}else if(Q->Inf==TRUE){
		EFp4_Copy(Ans, P);
		return;
	}else if(Fp4_Cmp(&P->x,&Q->x)==0&&Fp4_Cmp(&P->y,&Q->y)){
		EFp4_SetInf(Ans);
		return;
	}else if(Fp4_Cmp(&P->x,&Q->x)==0&&Fp4_Cmp(&P->y,&Q->y)==0){
		EFp4_ECD(Ans,P);
		return;
	}
	struct Fp4 lambda,tmp;
	struct EFp4 Ans_Copy;

	Fp4_Init(&lambda);
	Fp4_Init(&tmp);
	EFp4_Init(&Ans_Copy);

	Fp4_Sub(&tmp,&P->x,&Q->x);
	Fp4_Inv(&tmp,&tmp);
	Fp4_Sub(&lambda,&P->y,&Q->y);
	Fp4_Mul(&lambda,&lambda,&tmp);
	Fp4_Sub(&tmp,&P->x,&Q->x);

	Fp4_Sqr(&Ans_Copy.x,&lambda);
	Fp4_Sub(&Ans_Copy.x,&Ans_Copy.x,&P->x);
	Fp4_Sub(&Ans_Copy.x,&Ans_Copy.x,&Q->x);

	Fp4_Sub(&tmp,&P->x,&Ans_Copy.x);
	Fp4_Mul(&Ans_Copy.y,&lambda,&tmp);
	Fp4_Sub(&Ans_Copy.y,&Ans_Copy.y,&P->y);

	EFp4_Copy(Ans,&Ans_Copy);

	Fp4_Clear(&lambda);
	Fp4_Clear(&tmp);
	EFp4_Clear(&Ans_Copy);
}

void EFp4_ECD(struct EFp4 *Ans, struct EFp4 *P){
	if(P->Inf==TRUE){
		EFp4_SetInf(Ans);
		return;
	}else if(mpz_cmp_ui(P->y.tau[0],0)==0&&mpz_cmp_ui(P->y.tau[1],0)==0&&mpz_cmp_ui(P->y.tau[2],0)==0&&mpz_cmp_ui(P->y.tau[3],0)==0){
		EFp4_SetInf(Ans);
		return;
	}
	struct Fp4 lambda,tmp;
	struct EFp4 Ans_Copy;

	Fp4_Init(&lambda);
	Fp4_Init(&tmp);
	EFp4_Init(&Ans_Copy);

	Fp4_Mul(&lambda,&P->x,&P->x);
	Fp4_Add(&tmp,&lambda,&lambda);
	Fp4_Add(&lambda,&lambda,&tmp);
	Fp4_Add(&tmp,&P->y,&P->y);
	Fp4_Inv(&tmp,&tmp);
	Fp4_Mul(&lambda,&lambda,&tmp);

	Fp4_Mul(&Ans_Copy.x,&lambda,&lambda);
	Fp4_Add(&tmp,&P->x,&P->x);
	Fp4_Sub(&Ans_Copy.x,&Ans_Copy.x,&tmp);

	Fp4_Sub(&tmp,&P->x,&Ans_Copy.x);
	Fp4_Mul(&Ans_Copy.y,&lambda,&tmp);
	Fp4_Sub(&Ans_Copy.y,&Ans_Copy.y,&P->y);

	EFp4_Copy(Ans,&Ans_Copy);

	Fp4_Clear(&lambda);
	Fp4_Clear(&tmp);
	EFp4_Clear(&Ans_Copy);
}

void EFp4_SCM(struct EFp4 *Ans, struct EFp4 *P, mpz_t j){
	int i;
	int r;
	r = (int)mpz_sizeinbase(j,2);

	struct EFp4 Ans_Copy;
	EFp4_Init(&Ans_Copy);
	EFp4_Copy(&Ans_Copy,P);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp4_ECD(&Ans_Copy,&Ans_Copy);
			EFp4_ECA(&Ans_Copy,&Ans_Copy,P);
		}else{
			EFp4_ECD(&Ans_Copy,&Ans_Copy);
		}
	}

	EFp4_Copy(Ans,&Ans_Copy);
	EFp4_Clear(&Ans_Copy);
}

void EFp4_SCM2(struct EFp4 *Ans, struct EFp4 *P, mpz_t j){
	int TabSize = 2<<(WSize-1);
	struct EFp4 Table[TabSize],Ans_Copy;
	int i,r,k,Num;
	for(i=0;i<TabSize;i++){
		EFp4_Init(&Table[i]);
	}
	EFp4_Init(&Ans_Copy);
	
	EFp4_SetInf(&Ans_Copy);
	EFp4_SetInf(&Table[0]);
	for(i=1;i<TabSize;i++){
		EFp4_ECA(&Table[i],&Table[i-1],P);
	}
	
	//Size addjustment
	r=(int)mpz_sizeinbase(j,2);
	r+=(WSize-r%WSize)%WSize;
	r--;
	
	Num=0;
	for(k=1;k<WSize;k++){
		Num+=mpz_tstbit(j,r);
		Num=Num<<1;
		r--;
	}
	Num+=mpz_tstbit(j,r);
	r--;
	EFp4_ECA(&Ans_Copy,&Ans_Copy,&Table[Num]);
	for(i=r;i>=0;i-=WSize){
		for(k=0;k<WSize;k++){
			EFp4_ECD(&Ans_Copy,&Ans_Copy);
		}
		Num=0;
		for(k=1;k<WSize;k++){
			Num+=mpz_tstbit(j,r);
			Num=Num<<1;
			r--;
		}
		Num+=mpz_tstbit(j,r);
		r--;
		EFp4_ECA(&Ans_Copy,&Ans_Copy,&Table[Num]);
	}
	EFp4_Copy(Ans,&Ans_Copy);
	
	for(i=0;i<TabSize;i++){
		EFp4_Clear(&Table[i]);
	}
	EFp4_Clear(&Ans_Copy);
}

//-------------EFp8 functions--------------------------------
void EFp8_Init(struct EFp8 *A){
	Fp8_Init(&A->x);
	Fp8_Init(&A->y);
	A->Inf=FALSE;
}

void EFp8_Copy(struct EFp8 *A,struct EFp8 *B){
	Fp8_Copy(&A->x,&B->x);
	Fp8_Copy(&A->y,&B->y);
	A->Inf=B->Inf;
}

void EFp8_SetInf(struct EFp8 *A){
	Fp8_Init(&A->x);
	Fp8_Init(&A->y);
	A->Inf=TRUE;
}

void EFp8_Clear(struct EFp8 *A){
	Fp8_Clear(&A->x);
	Fp8_Clear(&A->y);
}

void EFp8_Rand(struct EFp8 *A){
	struct EFp8 Ans_Copy;
	struct Fp8 tmp,ECC_b;
	mpz_t j;
	EFp8_Init(&Ans_Copy);
	Fp8_Init(&tmp);
	Fp8_Init(&ECC_b);
	mpz_init(j);
	
	mpz_set_ui(j,3);
	Fp8_Set_ui(&ECC_b,ECCPara_b);
	while(Fp8_Legendre(&tmp)!=1){
		Fp8_Rand(&Ans_Copy.x);
		Fp8_Pow(&tmp,&Ans_Copy.x,j);
		Fp8_Add(&tmp,&tmp,&ECC_b);
	}
	Fp8_Sqrt(&Ans_Copy.y,&tmp);
	EFp8_Copy(A,&Ans_Copy);

	EFp8_Clear(&Ans_Copy);
	Fp8_Clear(&tmp);
	Fp8_Clear(&ECC_b);
	mpz_clear(j);
}

void EFp8_Show(struct EFp8 *A){
	gmp_printf("((%Zd,%Zd,%Zd,%Zd),(%Zd,%Zd,%Zd,%Zd))\n",A->x.theta[0].tau[0],A->x.theta[0].tau[1],A->x.theta[0].tau[2],A->x.theta[0].tau[3],A->x.theta[1].tau[0],A->x.theta[1].tau[1],A->x.theta[1].tau[2],A->x.theta[1].tau[3]);
	gmp_printf("((%Zd,%Zd,%Zd,%Zd),(%Zd,%Zd,%Zd,%Zd))\n",A->y.theta[0].tau[0],A->y.theta[0].tau[1],A->y.theta[0].tau[2],A->y.theta[0].tau[3],A->y.theta[1].tau[0],A->y.theta[1].tau[1],A->y.theta[1].tau[2],A->y.theta[1].tau[3]);
}

void EFp8_ECA(struct EFp8 *Ans, struct EFp8 *P, struct EFp8 *Q){
	if(P->Inf==TRUE){
		EFp8_Copy(Ans, Q);
		return;
	}else if(Q->Inf==TRUE){
		EFp8_Copy(Ans, P);
		return;
	}else if(Fp8_Cmp(&P->x,&Q->x)==0&&Fp8_Cmp(&P->y,&Q->y)){
		EFp8_SetInf(Ans);
		return;
	}else if(Fp8_Cmp(&P->x,&Q->x)==0&&Fp8_Cmp(&P->y,&Q->y)==0){
		EFp8_ECD(Ans,P);
		return;
	}

	struct Fp8 lambda,tmp;
	struct EFp8 Ans_Copy;

	Fp8_Init(&lambda);
	Fp8_Init(&tmp);
	EFp8_Init(&Ans_Copy);

	Fp8_Sub(&tmp,&P->x,&Q->x);
	Fp8_Inv(&tmp,&tmp);
	Fp8_Sub(&lambda,&P->y,&Q->y);
	Fp8_Mul(&lambda,&lambda,&tmp);

	Fp8_Mul(&Ans_Copy.x,&lambda,&lambda);
	Fp8_Sub(&Ans_Copy.x,&Ans_Copy.x,&P->x);
	Fp8_Sub(&Ans_Copy.x,&Ans_Copy.x,&Q->x);

	Fp8_Sub(&tmp,&P->x,&Ans_Copy.x);
	Fp8_Mul(&Ans_Copy.y,&lambda,&tmp);
	Fp8_Sub(&Ans_Copy.y,&Ans_Copy.y,&P->y);

	EFp8_Copy(Ans,&Ans_Copy);

	Fp8_Clear(&lambda);
	Fp8_Clear(&tmp);
	EFp8_Clear(&Ans_Copy);
}

void EFp8_ECD(struct EFp8 *Ans, struct EFp8 *P){
	if(P->Inf==TRUE){
		EFp8_SetInf(Ans);
		return;
	}else if(mpz_cmp_ui(P->y.theta[0].tau[0],0)==0&&mpz_cmp_ui(P->y.theta[0].tau[1],0)==0&&mpz_cmp_ui(P->y.theta[0].tau[2],0)==0&&mpz_cmp_ui(P->y.theta[0].tau[3],0)==0&&mpz_cmp_ui(P->y.theta[1].tau[0],0)==0&&mpz_cmp_ui(P->y.theta[1].tau[1],0)==0&&mpz_cmp_ui(P->y.theta[1].tau[2],0)==0&&mpz_cmp_ui(P->y.theta[1].tau[3],0)==0){
		EFp8_SetInf(Ans);
		return;
	}
	
	struct Fp8 lambda,tmp;
	struct EFp8 Ans_Copy;

	Fp8_Init(&lambda);
	Fp8_Init(&tmp);
	EFp8_Init(&Ans_Copy);

	Fp8_Mul(&lambda,&P->x,&P->x);
	Fp8_Add(&tmp,&lambda,&lambda);
	Fp8_Add(&lambda,&lambda,&tmp);
	Fp8_Add(&tmp,&P->y,&P->y);
	Fp8_Inv(&tmp,&tmp);
	Fp8_Mul(&lambda,&lambda,&tmp);

	Fp8_Mul(&Ans_Copy.x,&lambda,&lambda);
	Fp8_Add(&tmp,&P->x,&P->x);
	Fp8_Sub(&Ans_Copy.x,&Ans_Copy.x,&tmp);

	Fp8_Sub(&tmp,&P->x,&Ans_Copy.x);
	Fp8_Mul(&Ans_Copy.y,&lambda,&tmp);
	Fp8_Sub(&Ans_Copy.y,&Ans_Copy.y,&P->y);

	EFp8_Copy(Ans,&Ans_Copy);

	Fp8_Clear(&lambda);
	Fp8_Clear(&tmp);
	EFp8_Clear(&Ans_Copy);
}

void EFp8_SCM(struct EFp8 *Ans, struct EFp8 *P, mpz_t j){
	int i;
	int r;
	r = (int)mpz_sizeinbase(j,2);

	struct EFp8 Ans_Copy;
	EFp8_Init(&Ans_Copy);
	EFp8_Copy(&Ans_Copy,P);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp8_ECD(&Ans_Copy,&Ans_Copy);
			EFp8_ECA(&Ans_Copy,&Ans_Copy,P);
		}else{
			EFp8_ECD(&Ans_Copy,&Ans_Copy);
		}
	}

	EFp8_Copy(Ans,&Ans_Copy);
	EFp8_Clear(&Ans_Copy);
}

void EFp8_SCM2(struct EFp8 *Ans, struct EFp8 *P, mpz_t j){
	int TabSize = 2<<(WSize-1);
	struct EFp8 Table[TabSize],Ans_Copy;
	int i,r,k,Num;
	for(i=0;i<TabSize;i++){
		EFp8_Init(&Table[i]);
	}
	EFp8_Init(&Ans_Copy);
	
	EFp8_SetInf(&Ans_Copy);
	EFp8_SetInf(&Table[0]);
	for(i=1;i<TabSize;i++){
		EFp8_ECA(&Table[i],&Table[i-1],P);
	}
	
	//Size addjustment
	r=(int)mpz_sizeinbase(j,2);
	r+=(WSize-r%WSize)%WSize;
	r--;
	
	Num=0;
	for(k=1;k<WSize;k++){
		Num+=mpz_tstbit(j,r);
		Num=Num<<1;
		r--;
	}
	Num+=mpz_tstbit(j,r);
	r--;
	EFp8_ECA(&Ans_Copy,&Ans_Copy,&Table[Num]);
	for(i=r;i>=0;i-=WSize){
		for(k=0;k<WSize;k++){
			EFp8_ECD(&Ans_Copy,&Ans_Copy);
		}
		Num=0;
		for(k=1;k<WSize;k++){
			Num+=mpz_tstbit(j,r);
			Num=Num<<1;
			r--;
		}
		Num+=mpz_tstbit(j,r);
		r--;
		EFp8_ECA(&Ans_Copy,&Ans_Copy,&Table[Num]);
	}
	EFp8_Copy(Ans,&Ans_Copy);
	
	for(i=0;i<TabSize;i++){
		EFp8_Clear(&Table[i]);
	}
	EFp8_Clear(&Ans_Copy);
}

//-------------EFp24 functions--------------------------------
void EFp24_Init(struct EFp24 *A){
	Fp24_Init(&A->x);
	Fp24_Init(&A->y);
	A->Inf=FALSE;
}

void EFp24_Copy(struct EFp24 *A,struct EFp24 *B){
	Fp24_Copy(&A->x,&B->x);
	Fp24_Copy(&A->y,&B->y);
	A->Inf=B->Inf;
}

void EFp24_SetInf(struct EFp24 *A){
	Fp24_Init(&A->x);
	Fp24_Init(&A->y);
	A->Inf=TRUE;
}

void EFp24_Clear(struct EFp24 *A){
	Fp24_Clear(&A->x);
	Fp24_Clear(&A->y);
}

void EFp24_Rand(struct EFp24 *A){
	struct EFp24 Ans_Copy;
	struct Fp24 tmp,ECC_b;
	mpz_t j;
	EFp24_Init(&Ans_Copy);
	Fp24_Init(&tmp);
	Fp24_Init(&ECC_b);
	mpz_init(j);
	
	Fp24_Set_ui(&ECC_b,ECCPara_b);
	mpz_set_ui(j,3);
	while(Fp24_Legendre(&tmp)!=1){
		Fp24_Rand(&Ans_Copy.x);
		Fp24_Pow(&tmp,&Ans_Copy.x,j);
		Fp24_Add(&tmp,&tmp,&ECC_b);
	}
	Fp24_Sqrt(&Ans_Copy.y,&tmp);
	EFp24_Copy(A,&Ans_Copy);

	EFp24_Clear(&Ans_Copy);
	Fp24_Clear(&tmp);
	Fp24_Clear(&ECC_b);
	mpz_clear(j);
}

void EFp24_Show(struct EFp24 *A){
	printf("x=\n");
	gmp_printf("((%Zd,%Zd,%Zd,%Zd),(%Zd,%Zd,%Zd,%Zd))\n",A->x.omega[0].theta[0].tau[0],A->x.omega[0].theta[0].tau[1],A->x.omega[0].theta[0].tau[2],A->x.omega[0].theta[0].tau[3],A->x.omega[0].theta[1].tau[0],A->x.omega[0].theta[1].tau[1],A->x.omega[0].theta[1].tau[2],A->x.omega[0].theta[1].tau[3]);
	gmp_printf("((%Zd,%Zd,%Zd,%Zd),(%Zd,%Zd,%Zd,%Zd))\n",A->x.omega[1].theta[0].tau[0],A->x.omega[1].theta[0].tau[1],A->x.omega[1].theta[0].tau[2],A->x.omega[1].theta[0].tau[3],A->x.omega[1].theta[1].tau[0],A->x.omega[1].theta[1].tau[1],A->x.omega[1].theta[1].tau[2],A->x.omega[1].theta[1].tau[3]);
	gmp_printf("((%Zd,%Zd,%Zd,%Zd),(%Zd,%Zd,%Zd,%Zd))\n",A->x.omega[2].theta[0].tau[0],A->x.omega[2].theta[0].tau[1],A->x.omega[2].theta[0].tau[2],A->x.omega[2].theta[0].tau[3],A->x.omega[2].theta[1].tau[0],A->x.omega[2].theta[1].tau[1],A->x.omega[2].theta[1].tau[2],A->x.omega[2].theta[1].tau[3]);
	printf("y=\n");
	gmp_printf("((%Zd,%Zd,%Zd,%Zd),(%Zd,%Zd,%Zd,%Zd))\n",A->y.omega[0].theta[0].tau[0],A->y.omega[0].theta[0].tau[1],A->y.omega[0].theta[0].tau[2],A->y.omega[0].theta[0].tau[3],A->y.omega[0].theta[1].tau[0],A->y.omega[0].theta[1].tau[1],A->y.omega[0].theta[1].tau[2],A->y.omega[0].theta[1].tau[3]);
	gmp_printf("((%Zd,%Zd,%Zd,%Zd),(%Zd,%Zd,%Zd,%Zd))\n",A->y.omega[1].theta[0].tau[0],A->y.omega[1].theta[0].tau[1],A->y.omega[1].theta[0].tau[2],A->y.omega[1].theta[0].tau[3],A->y.omega[1].theta[1].tau[0],A->y.omega[1].theta[1].tau[1],A->y.omega[1].theta[1].tau[2],A->y.omega[1].theta[1].tau[3]);
	gmp_printf("((%Zd,%Zd,%Zd,%Zd),(%Zd,%Zd,%Zd,%Zd))\n",A->y.omega[2].theta[0].tau[0],A->y.omega[2].theta[0].tau[1],A->y.omega[2].theta[0].tau[2],A->y.omega[2].theta[0].tau[3],A->y.omega[2].theta[1].tau[0],A->y.omega[2].theta[1].tau[1],A->y.omega[2].theta[1].tau[2],A->y.omega[2].theta[1].tau[3]);
}

void EFp24_ECA(struct EFp24 *Ans, struct EFp24 *P, struct EFp24 *Q){
	if(P->Inf==TRUE){
		EFp24_Copy(Ans, Q);
		return;
	}else if(Q->Inf==TRUE){
		EFp24_Copy(Ans, P);
		return;
	}else if(Fp24_Cmp(&P->x,&Q->x)==0&&Fp24_Cmp(&P->y,&Q->y)){
		EFp24_SetInf(Ans);
		return;
	}else if(Fp24_Cmp(&P->x,&Q->x)==0&&Fp24_Cmp(&P->y,&Q->y)==0){
		EFp24_ECD(Ans,P);
		return;
	}

	struct Fp24 lambda,tmp;
	struct EFp24 Ans_Copy;

	Fp24_Init(&lambda);
	Fp24_Init(&tmp);
	EFp24_Init(&Ans_Copy);

	Fp24_Sub(&tmp,&P->x,&Q->x);
	Fp24_Inv(&tmp,&tmp);
	Fp24_Sub(&lambda,&P->y,&Q->y);
	Fp24_Mul(&lambda,&lambda,&tmp);

	Fp24_Mul(&Ans_Copy.x,&lambda,&lambda);
	Fp24_Sub(&Ans_Copy.x,&Ans_Copy.x,&P->x);
	Fp24_Sub(&Ans_Copy.x,&Ans_Copy.x,&Q->x);

	Fp24_Sub(&tmp,&P->x,&Ans_Copy.x);
	Fp24_Mul(&Ans_Copy.y,&lambda,&tmp);
	Fp24_Sub(&Ans_Copy.y,&Ans_Copy.y,&P->y);
	
	EFp24_Copy(Ans,&Ans_Copy);

	Fp24_Clear(&lambda);
	Fp24_Clear(&tmp);
	EFp24_Clear(&Ans_Copy);
}

void EFp24_ECD(struct EFp24 *Ans, struct EFp24 *P){
	mpz_t Zero;
	mpz_init(Zero);
	if(P->Inf==TRUE){
		EFp24_SetInf(Ans);
		mpz_clear(Zero);
		return;
	}else if(Fp24_Cmp_mpz(&P->y,Zero)==0){
		EFp24_SetInf(Ans);
		mpz_clear(Zero);
		return;
	}

	struct Fp24 lambda,tmp;
	struct EFp24 Ans_Copy;

	Fp24_Init(&lambda);
	Fp24_Init(&tmp);
	EFp24_Init(&Ans_Copy);

	Fp24_Mul(&lambda,&P->x,&P->x);
	Fp24_Add(&tmp,&lambda,&lambda);
	Fp24_Add(&lambda,&lambda,&tmp);
	Fp24_Add(&tmp,&P->y,&P->y);
	Fp24_Inv(&tmp,&tmp);
	Fp24_Mul(&lambda,&lambda,&tmp);

	Fp24_Mul(&Ans_Copy.x,&lambda,&lambda);
	Fp24_Add(&tmp,&P->x,&P->x);
	Fp24_Sub(&Ans_Copy.x,&Ans_Copy.x,&tmp);

	Fp24_Sub(&tmp,&P->x,&Ans_Copy.x);
	Fp24_Mul(&Ans_Copy.y,&lambda,&tmp);
	Fp24_Sub(&Ans_Copy.y,&Ans_Copy.y,&P->y);

	EFp24_Copy(Ans,&Ans_Copy);

	Fp24_Clear(&lambda);
	Fp24_Clear(&tmp);
	EFp24_Clear(&Ans_Copy);
}
/*
void EFp24_SCM(struct EFp24 *Ans, struct EFp24 *P, mpz_t j){
	int i;
	int r;
	r = (int)mpz_sizeinbase(j,2);

	struct EFp24 Ans_Copy;
	EFp24_Init(&Ans_Copy);
	EFp24_Copy(&Ans_Copy,P);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp24_ECD(&Ans_Copy,&Ans_Copy);
			EFp24_ECA(&Ans_Copy,&Ans_Copy,P);
		}else{
			EFp24_ECD(&Ans_Copy,&Ans_Copy);
		}
	}

	EFp24_Copy(Ans,&Ans_Copy);
	EFp24_Clear(&Ans_Copy);
}*/

//Window Method
void EFp24_SCM(struct EFp24 *Ans, struct EFp24 *P, mpz_t j){
	int TabSize = 2<<(WSize-1);
	struct EFp24 Table[TabSize],Ans_Copy;
	int i,r,k,Num;
	for(i=0;i<TabSize;i++){
		EFp24_Init(&Table[i]);
	}
	EFp24_Init(&Ans_Copy);
	
	EFp24_SetInf(&Ans_Copy);
	EFp24_SetInf(&Table[0]);
	for(i=1;i<TabSize;i++){
		EFp24_ECA(&Table[i],&Table[i-1],P);
	}
	
	//Size addjustment
	r=(int)mpz_sizeinbase(j,2);
	r+=(WSize-r%WSize)%WSize;
	r--;
	
	Num=0;
	for(k=1;k<WSize;k++){
		Num+=mpz_tstbit(j,r);
		Num=Num<<1;
		r--;
	}
	Num+=mpz_tstbit(j,r);
	r--;
	EFp24_ECA(&Ans_Copy,&Ans_Copy,&Table[Num]);
	for(i=r;i>=0;i-=WSize){
		for(k=0;k<WSize;k++){
			EFp24_ECD(&Ans_Copy,&Ans_Copy);
		}
		Num=0;
		for(k=1;k<WSize;k++){
			Num+=mpz_tstbit(j,r);
			Num=Num<<1;
			r--;
		}
		Num+=mpz_tstbit(j,r);
		r--;
		EFp24_ECA(&Ans_Copy,&Ans_Copy,&Table[Num]);
	}
	EFp24_Copy(Ans,&Ans_Copy);
	
	for(i=0;i<TabSize;i++){
		EFp24_Clear(&Table[i]);
	}
	EFp24_Clear(&Ans_Copy);
}

void EFp24_Inv(struct EFp24 *Ans, struct EFp24 *P){
	Fp24_Copy(&Ans->x,&P->x);
	Fp24_Neg(&Ans->y,&P->y);
}

void EFp24_Frob(struct EFp24 *Ans, struct EFp24 *P){
	Fp24_Frobp(&Ans->x,&P->x);
	Fp24_Frobp(&Ans->y,&P->y);
}
/*
void EFp24_G1_Rand(struct EFp24 *A){
	struct EFp24 Ans_Copy;
	struct EFp Rand;
	EFp24_Init(&Ans_Copy);
	EFp_Init(&Rand);

	EFp_Rand(&Rand);
	
	mpz_t i;
	mpz_init(i);
	mpz_set_ui(i,14700131460293907);
	EFp_SCM(&Rand,&Rand,i);//Rand is in {E[r] and E(Fp)}

	mpz_set(Ans_Copy.x.omega[0].theta[0].tau[0],Rand.x);
	mpz_set(Ans_Copy.x.omega[0].theta[0].tau[1],Rand.x);
	mpz_set(Ans_Copy.x.omega[0].theta[0].tau[2],Rand.x);
	mpz_set(Ans_Copy.x.omega[0].theta[0].tau[3],Rand.x);
	mpz_set(Ans_Copy.x.omega[1].theta[0].tau[0],Rand.x);
	mpz_set(Ans_Copy.x.omega[1].theta[0].tau[1],Rand.x);
	mpz_set(Ans_Copy.x.omega[1].theta[0].tau[2],Rand.x);
	mpz_set(Ans_Copy.x.omega[1].theta[0].tau[3],Rand.x);
	mpz_set(Ans_Copy.x.omega[2].theta[0].tau[0],Rand.x);
	mpz_set(Ans_Copy.x.omega[2].theta[0].tau[1],Rand.x);
	mpz_set(Ans_Copy.x.omega[2].theta[0].tau[2],Rand.x);
	mpz_set(Ans_Copy.x.omega[2].theta[0].tau[3],Rand.x);
	mpz_set(Ans_Copy.y.omega[0].theta[0].tau[0],Rand.y);
	mpz_set(Ans_Copy.y.omega[0].theta[0].tau[1],Rand.y);
	mpz_set(Ans_Copy.y.omega[0].theta[0].tau[2],Rand.y);
	mpz_set(Ans_Copy.y.omega[0].theta[0].tau[3],Rand.y);
	mpz_set(Ans_Copy.y.omega[1].theta[0].tau[0],Rand.y);
	mpz_set(Ans_Copy.y.omega[1].theta[0].tau[1],Rand.y);
	mpz_set(Ans_Copy.y.omega[1].theta[0].tau[2],Rand.y);
	mpz_set(Ans_Copy.y.omega[1].theta[0].tau[3],Rand.y);
	mpz_set(Ans_Copy.y.omega[2].theta[0].tau[0],Rand.y);
	mpz_set(Ans_Copy.y.omega[2].theta[0].tau[1],Rand.y);
	mpz_set(Ans_Copy.y.omega[2].theta[0].tau[2],Rand.y);
	mpz_set(Ans_Copy.y.omega[2].theta[0].tau[3],Rand.y);

	EFp24_Copy(A,&Ans_Copy);
	EFp24_Clear(&Ans_Copy);
	EFp_Clear(&Rand);
	mpz_clear(i);
}
*/
void EFp24_G1_Rand(struct EFp24 *A){
	struct EFp24 Ans_Copy,P_Frob;
	struct Fp24 tmp,ECC_b;
	mpz_t j,R;
	EFp24_Init(&Ans_Copy);
	EFp24_Init(&P_Frob);
	Fp24_Init(&tmp);
	Fp24_Init(&ECC_b);
	mpz_init(j);
	mpz_init(R);
	
	mpz_tdiv_q(R,r24,order);
	mpz_tdiv_q(R,R,order);

	Fp24_Set_ui(&ECC_b,ECCPara_b);
	mpz_set_ui(j,3);
	while(Fp24_Legendre(&tmp)!=1){
		Fp24_Rand(&Ans_Copy.x);
		Fp24_Pow(&tmp,&Ans_Copy.x,j);
		Fp24_Add(&tmp,&tmp,&ECC_b);
	}
	Fp24_Sqrt(&Ans_Copy.y,&tmp);
	EFp24_SCM(&Ans_Copy,&Ans_Copy,R);

	Fp24_Frobp(&P_Frob.x,&Ans_Copy.x);
	Fp24_Frobp(&P_Frob.y,&Ans_Copy.y);
	EFp24_SCM(&Ans_Copy,&Ans_Copy,prime);
	Fp24_Neg(&Ans_Copy.y,&Ans_Copy.y);
	EFp24_ECA(&Ans_Copy,&P_Frob,&Ans_Copy);//Frob(R)-[p]R
	EFp24_Copy(A,&Ans_Copy);

	EFp24_Clear(&Ans_Copy);
	Fp24_Clear(&tmp);
	Fp24_Clear(&ECC_b);
	mpz_clear(j);
	mpz_clear(R);
}

void EFp24_OptG1_Rand(struct EFp24 *A){
	struct EFp24 Ans_Copy;
	struct EFp8 Rand;
	EFp24_Init(&Ans_Copy);
	EFp8_Init(&Rand);

	EFp8_Rand(&Rand);
	
	Fp8_Copy(&Ans_Copy.x.omega[0],&Rand.x);
	Fp8_Copy(&Ans_Copy.x.omega[1],&Rand.x);
	Fp8_Copy(&Ans_Copy.x.omega[2],&Rand.x);
	Fp8_Copy(&Ans_Copy.y.omega[0],&Rand.y);
	Fp8_Copy(&Ans_Copy.y.omega[1],&Rand.y);
	Fp8_Copy(&Ans_Copy.y.omega[2],&Rand.y);
	EFp24_Copy(A,&Ans_Copy);
	EFp24_Clear(&Ans_Copy);
	EFp8_Clear(&Rand);
}

void EFp24_G2_Rand(struct EFp24 *Ans){
	struct EFp24 Ans_Copy,P,P_Frob;
	mpz_t R;
	EFp24_Init(&Ans_Copy);
	EFp24_Init(&P);
	EFp24_Init(&P_Frob);
	mpz_init(R);

	mpz_tdiv_q(R,r24,order);
	mpz_tdiv_q(R,R,order);
	
	EFp24_Rand(&P);
	EFp24_SCM(&P,&P,R);
	Fp24_Frobp(&P_Frob.x,&P.x);
	Fp24_Frobp(&P_Frob.y,&P.y);
	Fp24_Neg(&Ans_Copy.y,&P.y);
	Fp24_Copy(&Ans_Copy.x,&P.x);
	//EFp24_SCM(&Ans_Copy,&Ans_Copy,prime);

	EFp24_ECA(&Ans_Copy,&Ans_Copy,&P_Frob);
	EFp24_Copy(Ans,&Ans_Copy);
	
	EFp24_Clear(&Ans_Copy);
	EFp24_Clear(&P);
	EFp24_Clear(&P_Frob);
	mpz_clear(R);
}
/*
void EFp24_Check(){
	struct EFp24 A,B;
	EFp24_Init(&A);
	EFp24_Init(&B);
	mpz_t j;
	mpz_init(j);
	
	EFp24_Rand(&A);
	EFp24_Show(&A);
	mpz_sub_ui(j,prime,49);
	EFp24_SCM(&B,&A,j);
	EFp24_Show(&B);
}*/

//-------------Miller Algorithm--------------------------------
void Miller_Algo(struct Fp24 *Ans, struct EFp24 *P, struct EFp24 *Q, mpz_t roop){
	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(roop,2);

	struct Fp24 lambda,l,l_sum,x,y,v,v_sum,tmp;
	struct EFp24 T;
	Fp24_Init(&lambda);
	Fp24_Init(&l);
	Fp24_Init(&l_sum);
	Fp24_Init(&x);
	Fp24_Init(&y);
	Fp24_Init(&v);
	Fp24_Init(&v_sum);
	Fp24_Init(&tmp);
	EFp24_Init(&T);
	
	EFp24_Copy(&T,P);
	Fp24_Set1(&l_sum);
	Fp24_Set1(&v_sum);

	for(i=r-2;i>=0;i--){
		Fp24_Sqr(&l_sum,&l_sum);
		Fp24_Sqr(&v_sum,&v_sum);
		
		Fp24_Sqr(&lambda,&T.x);
		Fp24_Add(&tmp,&lambda,&lambda);
		Fp24_Add(&lambda,&tmp,&lambda);
		Fp24_Add(&tmp,&T.y,&T.y);
		Fp24_Inv(&tmp,&tmp);
		Fp24_Mul(&lambda,&lambda,&tmp);	//lambda <- (3x_T^2+a)/2y_T
		
		Fp24_Sub(&l,&Q->x,&T.x);
		Fp24_Mul(&l,&l,&lambda);
		Fp24_Sub(&tmp,&Q->y,&T.y);
		Fp24_Sub(&l,&l,&tmp);	//l <- (x_Q-x_T)*lambda-(y_Q-y_T)
	
		Fp24_Sqr(&x,&lambda);
		Fp24_Add(&tmp,&T.x,&T.x);
		Fp24_Sub(&x,&x,&tmp);	//x <- lambda^2-2x_T

		Fp24_Sub(&y,&T.x,&x);
		Fp24_Mul(&y,&y,&lambda);
		Fp24_Sub(&y,&y,&T.y);	//y <- (x_T-x)*lambda-y_T
	
		Fp24_Sub(&v,&Q->x,&x);	//v <- x_Q-x

		Fp24_Copy(&T.x,&x);	//R <- 2T
		Fp24_Copy(&T.y,&y);

		Fp24_Mul(&l_sum,&l_sum,&l);
		Fp24_Mul(&v_sum,&v_sum,&v);

		if(mpz_tstbit(roop,i)==1){
			Fp24_Sub(&lambda,&P->y,&T.y);
			Fp24_Sub(&tmp,&P->x,&T.x);
			Fp24_Inv(&tmp,&tmp);
			Fp24_Mul(&lambda,&lambda,&tmp);		//lambda <- (y_P-y_T)/(x_P-x_T)

			Fp24_Sub(&l,&Q->x,&P->x);
			Fp24_Mul(&l,&l,&lambda);
			Fp24_Sub(&tmp,&Q->y,&P->y);
			Fp24_Sub(&l,&l,&tmp);	//l <- (x_Q-x_P)*lambda-(y_Q-y_P)
	
			Fp24_Sqr(&x,&lambda);
			Fp24_Sub(&x,&x,&T.x);
			Fp24_Sub(&x,&x,&P->x);	//x <- lambda^2-x_t-x_P

			Fp24_Sub(&y,&P->x,&x);
			Fp24_Mul(&y,&y,&lambda);
			Fp24_Sub(&y,&y,&P->y);	//y <- (x_P-x)*lambda-y_P

			Fp24_Sub(&v,&Q->x,&x);	//v <- x_Q-x

			Fp24_Copy(&T.x,&x);		//R <- (T+P)=(x,y)
			Fp24_Copy(&T.y,&y);

			Fp24_Mul(&l_sum,&l_sum,&l);
			Fp24_Mul(&v_sum,&v_sum,&v);
		}
	}

	Fp24_Inv(&v_sum,&v_sum);
	Fp24_Mul(Ans,&l_sum,&v_sum);
	
	Fp24_Clear(&lambda);
	Fp24_Clear(&l);
	Fp24_Clear(&l_sum);
	Fp24_Clear(&x);
	Fp24_Clear(&y);
	Fp24_Clear(&v);
	Fp24_Clear(&v_sum);
	Fp24_Clear(&tmp);
	EFp24_Clear(&T);
}
/*
void Miller_Algo(struct Fp24 *Ans, struct EFp24 *P, struct EFp24 *Q, mpz_t roop){
	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(roop,2);

	struct Fp24 f,tmp;
	struct EFp24 T;
	Fp24_Init(&f);
	Fp24_Init(&tmp);
	EFp24_Init(&T);
	
	Fp24_Set1(&f);
	EFp24_Copy(&T,P);

	struct EFp24 X;
	EFp24_Init(&X);

	for(i=r-2;i>=0;i--){
		Fp24_Sqr(&f,&f);
		FDBL(&tmp,&T,&T,Q);
		Fp24_Mul(&f,&f,&tmp);
		if(mpz_tstbit(roop,i)==1){
			FADD(&tmp,&T,&T,P,Q);
			Fp24_Mul(&f,&f,&tmp);
		}
	}
	Fp24_Copy(Ans,&f);
	
	Fp24_Clear(&f);
	Fp24_Clear(&tmp);
	EFp24_Clear(&T);
}
*/
void FDBL(struct Fp24 *Ans, struct EFp24 *R, struct EFp24 *T, struct EFp24 *Q){
	struct Fp24 lambda,l,x,y,v,g,tmp;
	Fp24_Init(&lambda);
	Fp24_Init(&l);
	Fp24_Init(&x);
	Fp24_Init(&y);
	Fp24_Init(&v);
	Fp24_Init(&g);
	Fp24_Init(&tmp);
	
	Fp24_Sqr(&lambda,&T->x);
	Fp24_Add(&tmp,&lambda,&lambda);
	Fp24_Add(&lambda,&tmp,&lambda);
	Fp24_Add(&tmp,&T->y,&T->y);
	Fp24_Inv(&tmp,&tmp);
	Fp24_Mul(&lambda,&lambda,&tmp);	//lambda <- (3x_T^2+a)/2y_T

	Fp24_Sub(&l,&Q->x,&T->x);
	Fp24_Mul(&l,&l,&lambda);
	Fp24_Sub(&tmp,&Q->y,&T->y);
	Fp24_Sub(&l,&l,&tmp);	//l <- (x_Q-x_T)*lambda-(y_Q-y_T)

	Fp24_Sqr(&x,&lambda);
	Fp24_Add(&tmp,&T->x,&T->x);
	Fp24_Sub(&x,&x,&tmp);	//x <- lambda^2-2x_T

	Fp24_Sub(&y,&T->x,&x);
	Fp24_Mul(&y,&y,&lambda);
	Fp24_Sub(&y,&y,&T->y);	//y <- (x_T-x)*lambda-y_T
	
	Fp24_Sub(&v,&Q->x,&x);	//v <- x_Q-x
	Fp24_Inv(&tmp,&v);
	Fp24_Mul(&g,&l,&tmp);	//g <- l*v^(-1)

	Fp24_Copy(&R->x,&x);	//R <- 2T
	Fp24_Copy(&R->y,&y);
	Fp24_Copy(Ans,&g);


	//Fp24_Show(&l);
	//printf("\n");

	Fp24_Clear(&lambda);
	Fp24_Clear(&l);
	Fp24_Clear(&x);
	Fp24_Clear(&y);
	Fp24_Clear(&v);
	Fp24_Clear(&g);
	Fp24_Clear(&tmp);
}

void FADD(struct Fp24 *Ans, struct EFp24 *R, struct EFp24 *T, struct EFp24 *P, struct EFp24 *Q){
	struct Fp24 lambda,l,x,y,v,g,tmp;
	Fp24_Init(&lambda);
	Fp24_Init(&l);
	Fp24_Init(&x);
	Fp24_Init(&y);
	Fp24_Init(&v);
	Fp24_Init(&g);
	Fp24_Init(&tmp);
	
	Fp24_Sub(&lambda,&P->y,&T->y);
	Fp24_Sub(&tmp,&P->x,&T->x);
	Fp24_Inv(&tmp,&tmp);
	Fp24_Mul(&lambda,&lambda,&tmp);		//lambda <- (y_P-y_T)/(x_P-x_T)

	Fp24_Sub(&l,&Q->x,&P->x);
	Fp24_Mul(&l,&l,&lambda);
	Fp24_Sub(&tmp,&Q->y,&P->y);
	Fp24_Sub(&l,&l,&tmp);	//l <- (x_Q-x_P)*lambda-(y_Q-y_P)
	
	Fp24_Sqr(&x,&lambda);
	Fp24_Sub(&x,&x,&T->x);
	Fp24_Sub(&x,&x,&P->x);	//x <- lambda^2-x_t-x_P

	Fp24_Sub(&y,&P->x,&x);
	Fp24_Mul(&y,&y,&lambda);
	Fp24_Sub(&y,&y,&P->y);	//y <- (x_P-x)*lambda-y_P

	Fp24_Sub(&v,&Q->x,&x);
	Fp24_Inv(&tmp,&v);
	Fp24_Mul(&g,&l,&tmp);	//g <- l*v^(-1)

	Fp24_Copy(&R->x,&x);		//R <- 2T=(x,y)
	Fp24_Copy(&R->y,&y);
	Fp24_Copy(Ans,&g);

	//Fp24_Show(&l);
	//printf("\n");
	Fp24_Clear(&lambda);
	Fp24_Clear(&l);
	Fp24_Clear(&x);
	Fp24_Clear(&y);
	Fp24_Clear(&v);
	Fp24_Clear(&g);
	Fp24_Clear(&tmp);
}

void Final_Exp(struct Fp24 *Ans,struct Fp24 *A){
	mpz_t tmp;
	mpz_init(tmp);

	mpz_pow_ui(tmp,prime,24);
	mpz_sub_ui(tmp,tmp,1);
	mpz_div(tmp,tmp,order);

	Fp24_Pow(Ans,A,tmp);

	mpz_clear(tmp);
}

void Tate_Pairing(struct Fp24 *Ans, struct EFp24 *G1, struct EFp24 *G2){
	struct Fp24 Ans_Copy;
	Fp24_Init(&Ans_Copy);

	Miller_Algo(&Ans_Copy,G1,G2,order);
	printf("Mill\n");
	Fp24_Show(&Ans_Copy);
	Final_Exp(&Ans_Copy,&Ans_Copy);
	printf("FinalExp\n");
	Fp24_Show(&Ans_Copy);
	Fp24_Copy(Ans,&Ans_Copy);

	Fp24_Clear(&Ans_Copy);
}

void Ate_Pairing(struct Fp24 *Ans, struct EFp24 *G1, struct EFp24 *G2){
	struct Fp24 Ans_Copy;
	Fp24_Init(&Ans_Copy);

	Miller_Algo(&Ans_Copy,G1,G2,z);
	Final_Exp(&Ans_Copy,&Ans_Copy);
	Fp24_Copy(Ans,&Ans_Copy);

	Fp24_Clear(&Ans_Copy);
}

void Mill_SCM_Check(){
	struct EFp24 A,B;
	EFp24_Init(&A);
	EFp24_Init(&B);
	struct Fp24 D;
	Fp24_Init(&D);

	EFp24_Rand(&A);
	EFp24_G2_Rand(&B);

	Tate_Pairing(&D,&A,&B);
	EFp24_SCM(&A,&A,z);
	EFp24_Show(&A);

	EFp24_Clear(&A);
	EFp24_Clear(&B);
	Fp24_Clear(&D);
}

void Tate_Check(){
	struct EFp24 P,Q;
	struct Fp24 tmp;
	EFp24_Init(&P);
	EFp24_Init(&Q);
	Fp24_Init(&tmp);
	mpz_t a,b,c;
	mpz_init(a);
	mpz_init(b);
	mpz_init(c);
	mpz_set_ui(a,2);
	mpz_set_ui(b,3);
	mpz_set_ui(c,6);

	EFp24_G1_Rand(&P);
	EFp24_G2_Rand(&Q);
	printf("G1\n");
	EFp24_Show(&P);
	printf("G2\n");
	EFp24_Show(&Q);
	Tate_Pairing(&tmp,&P,&Q);
	Fp24_Pow(&tmp,&tmp,c);
	printf("result\n");
	Fp24_Show(&tmp);
	EFp24_SCM(&P,&P,a);
	EFp24_SCM(&Q,&Q,b);
	Tate_Pairing(&tmp,&P,&Q);
	printf("\n");
	printf("result\n");
	Fp24_Show(&tmp);
	
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(c);
	EFp24_Clear(&P);
	EFp24_Clear(&Q);
	Fp24_Clear(&tmp);
}

void Tate_Test(){
	clock_t start,end;
	struct EFp24 P,Q;
	struct Fp24 tmp;
	EFp24_Init(&P);
	EFp24_Init(&Q);
	Fp24_Init(&tmp);

	EFp24_G1_Rand(&P);
	EFp24_G2_Rand(&Q);
	printf("G1\n");
	EFp24_Show(&P);
	printf("G2\n");
	EFp24_Show(&Q);
	start = clock();
	Tate_Pairing(&tmp,&P,&Q);
	end = clock();
	printf("result\n");
	Fp24_Show(&tmp);//hoge

	printf("TatePairing Time %f\n",(double)(end-start)/CLOCKS_PER_SEC);

	EFp24_Clear(&P);
	EFp24_Clear(&Q);
	Fp24_Clear(&tmp);

}

void Ate_Check(){
	struct EFp24 P,Q;
	struct Fp24 tmp;
	EFp24_Init(&P);
	EFp24_Init(&Q);
	Fp24_Init(&tmp);
	mpz_t a,b,c;
	mpz_init(a);
	mpz_init(b);
	mpz_init(c);
	mpz_set_ui(a,2);
	mpz_set_ui(b,3);
	mpz_set_ui(c,6);

	EFp24_G1_Rand(&P);
	EFp24_G2_Rand(&Q);
	Ate_Pairing(&tmp,&Q,&P);
	Fp24_Pow(&tmp,&tmp,c);
	Fp24_Show(&tmp);
	EFp24_SCM(&P,&P,a);
	EFp24_SCM(&Q,&Q,b);
	Ate_Pairing(&tmp,&Q,&P);
	printf("\n");
	Fp24_Show(&tmp);
	
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(c);
	EFp24_Clear(&P);
	EFp24_Clear(&Q);
	Fp24_Clear(&tmp);
}

void Set_Parameter(mpz_t z,mpz_t prime,mpz_t order,mpz_t trace){
	mpz_t p4,t4,p8,t8,p24,t24,tmp1,tmp2;
	mpz_init(p4);
	mpz_init(t4);
	mpz_init(p8);
	mpz_init(t8);
	mpz_init(p24);
	mpz_init(t24);
	mpz_init(tmp1);
	mpz_init(tmp2);
	
	mpz_pow_ui(tmp1,z,8);
	mpz_add_ui(tmp1,tmp1,1);
	mpz_pow_ui(tmp2,z,4);
	mpz_sub(tmp1,tmp1,tmp2);
	mpz_set(order,tmp1);
	mpz_sub_ui(tmp2,z,1);
	mpz_mul(tmp2,tmp2,tmp2);
	mpz_mul(tmp1,tmp1,tmp2);
	mpz_div_ui(tmp1,tmp1,3);
	mpz_add(tmp1,tmp1,z);
	mpz_set(prime,tmp1);
	mpz_add_ui(trace,z,1);

	mpz_sub(r,prime,trace);
	mpz_add_ui(r,r,1);

	mpz_mul(t4,trace,trace);
	mpz_mul_ui(tmp2,prime,2);
	mpz_sub(t4,t4,tmp2);
	mpz_mul(t4,t4,t4);
	
	mpz_mul(tmp2,prime,prime);
	mpz_mul(p4,tmp2,tmp2);	//p4=p^4
	mpz_mul_ui(tmp2,tmp2,2);
	mpz_sub(t4,t4,tmp2);	//t4=(t^2-2p)^2-2p^2
	
	mpz_sub(r4,p4,t4);
	mpz_add_ui(r4,r4,1);

	mpz_mul(t8,t4,t4);
	mpz_mul_ui(tmp2,p4,2);
	mpz_sub(t8,t8,tmp2);
	mpz_mul(p8,p4,p4);

	mpz_sub(r8,p8,t8);
	mpz_add_ui(r8,r8,1);

	mpz_pow_ui(t24,t8,3);
	mpz_mul_ui(tmp2,p8,3);
	mpz_mul(tmp2,tmp2,t8);
	mpz_sub(t24,t24,tmp2);
	mpz_pow_ui(p24,p8,3);

	mpz_sub(r24,p24,t24);
	mpz_add_ui(r24,r24,1);	//r24 = order of Fp24

	ECCPara_b = 4;

	mpz_clear(p4);
	mpz_clear(t4);
	mpz_clear(p8);
	mpz_clear(t8);
	mpz_clear(p24);
	mpz_clear(t24);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
}

void ParaSerch(){
	mpz_t p,r,t,tmp,tmp2;
	mpz_init(p);
	mpz_init(r);
	mpz_init(t);
	mpz_init(tmp);
	mpz_init(tmp2);
	mpz_t a,b,c;
	mpz_init(a);
	mpz_init(b);
	mpz_init(c);
	mpz_set_str(z,"31010192",10);
	mpz_set_ui(a,0);

	mpz_t d,e,f,g;
	mpz_init(d);
	mpz_init(e);
	mpz_init(f);
	mpz_init(g);
	while(mpz_cmp_ui(a,1)!=0||mpz_cmp_ui(e,1)==0||mpz_cmp_ui(g,1)==0||mpz_cmp_ui(tmp,3)!=0){
	mpz_add_ui(z,z,1);
	Set_Parameter(z,p,r,t);
	mpz_mod_ui(c,p,7);
	gmp_printf("z=%Zd\np=%Zd\nr=%Zd\nt=%Zd\np%7=%Zd\n",z,p,r,t,c);
	printf("r=");
	mpz_out_str(stdout,2,r);
	printf("\n");
	mpz_set(prime,p);
	mpz_set_ui(a,3);
	Fp_Inv(b,a);
	mpz_mul(a,a,b);
	mpz_mod(a,a,prime);
	gmp_printf("a=%Zd\n",a);
	
	mpz_set_ui(d,61);
	mpz_sub_ui(e,prime,1);
	mpz_tdiv_q_ui(e,e,2);
	mpz_powm(e,d,e,prime);
	mpz_set_ui(f,61);
	mpz_sub_ui(g,prime,1);
	mpz_tdiv_q_ui(g,g,2);
	mpz_powm(g,f,g,prime);

	mpz_mod_ui(tmp,p,7);
	}
	
	mpz_sub_ui(tmp,p,1);
	mpz_div_ui(tmp,tmp,3);
	mpz_set_ui(tmp2,11);
	Fp_Pow(tmp2,tmp2,tmp);
	mpz_mod(tmp2,tmp2,p);
	gmp_printf("%Zd\n",tmp2);
	mpz_set(order,r);
	mpz_set(trace,t);
	//Tate_Check();
	mpz_clear(p);
	mpz_clear(r);
	mpz_clear(t);
	mpz_clear(tmp);
	mpz_clear(tmp2);

}

//-----------------------------------main-------------------------------------------
//bit length
//prime : 275
//p^k : 5935
//r : 222
int main(void){
	unsigned int now = (unsigned int)time(0);
	srand(now);

	//---------EFp Parameter set-----------
	//mpz_set_str(z,"31008931",10);
	mpz_set_str(z,"210000940",10);
	Set_Parameter(z,prime,order,trace);
	gmp_printf("z=%Zd\nprime=%Zd\norder=%Zd\ntrace=%Zd\n\n",z,prime,order,trace);
	//mpz_out_str(stdout,2,prime);
	//printf("\n");

	Tate_Test();
	//ParaSerch();
	//Tate_Check();
	//Ate_Check();

	
	return 0;
}
