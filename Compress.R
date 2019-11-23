Compressing Code in R
#library(ape);
library(stringr);
library(gtools)
#Program: DNA coding.
ThreeWordSeq = function(seq) 
{
x=c("A","C","G","T");
valores =permutations(n=4,r=3,v=x,repeats.allowed=T);

seq3=getNumSeq(seq);
tripletas=rep(0,64);
return(tripletas);
j=1;
while(j<=length(seq))
	{		
	char=paste0(seq[i],seq[i+1],seq[i+2]);
	char=toupper(chr);
	if(suma==350){valor=76}

	j=j+3;
	}
}
TwoWordSeq = function(seq) 
	{
	largo=length(seq);
	largo=largo-1;
	AA=0;CA=0;GA=0;TA=0;AC=0;CC=0;GC=0;TC=0;AG=0;CG=0;GG=0;TG=0;AT=0;CT=0;GT=0;TT=0;
	for (i in 1:largo) 
		{
		chr=paste0(seq[i],seq[i+1]);
		chr=toupper(chr);
			if (chr=="AA") {
				AA = AA+1;	
			}else if (chr=="AC") {
				AC=AC+1;	
			}else if (chr=="AG") {
				AG=AG+1;
			}else if (chr=="AT") {
				AT=AT+1;
			}else if (chr=="CA") {
				CA=CA+1;
			}else if (chr=="CC") {
				CC=CC+1;	
			}else if (chr=="CG") {
				CG=CG+1;
			}else if (chr=="CT") {
				CT=CT+1;			
			}else if (chr=="GA") {
				GA=GA+1;
			}else if (chr=="GC") {
				GC=GC+1;	
			}else if (chr=="GG") {
				GG=GG+1;
			}else if (chr=="GT") {
				GT=GT+1;	
			}else if (chr=="TA") {
				TA=TA+1;
			}else if (chr=="TC") {
				TC=TC+1;	
			}else if (chr=="TG") {
				TG=TG+1;
			}else if (chr=="TT") {
				TT=TT+1;	
			}
		}
	TwoWord=c(AA,AC,AG,AT);
	TwoWord=c(TwoWord,c(CA,CC,CG,CT));
	TwoWord=c(TwoWord,c(GA,GC,GG,GT));
	TwoWord=c(TwoWord,c(TA,TC,TG,TT));
	return(TwoWord);
	}
	
#simulated sequence using iid model
x= c("a","c","g","t");
seq2 = sample(x, 10000, replace=TRUE);
TwoWordMatrix = TwoWordSeq(seq2);

	#Function transfer DNA seq to numerical: a=1; c=2; g=3; t=4;
	#Input: a sequence (array) of "a", "c", "g", "t"
	#Output: the corresponding numeric values
	getNumSeq = function(seq) {
		numericalSeq = array(0, length(seq));
		for (i in 1:length(seq)) {
			if (seq[i]=="a") {
				numericalSeq[i] = 1;	
			}else if (seq[i]=="c") {
				numericalSeq[i] = 2;	
			}else if (seq[i]=="g") {
				numericalSeq[i] = 3;	
			}else if (seq[i]=="t") {
				numericalSeq[i] = 4;	
			}
		}
		return(numericalSeq);
	}
	#source();

is.even <- function(x) x %% 2 == 0
is.odd <- function(x) x %% 2 != 0

codificador = function(sequency) 
{
contador=0;
calibrador=rep(0,34);
sumatoria=rep(0,34);
codificacion=rep("0",21);
codificacion[18]="!";
codificacion[20]="!";

j=1;

while (j<=length(sequency)) 
	{
	contador=contador+1;
	sumatoria[contador]=strtoi(sequency[j])+strtoi(sequency[j+1])-strtoi(sequency[j+2]);
	calibrador[contador]=(2*(strtoi(sequency[j])*strtoi(sequency[j+2])));
      calibrador[contador] =calibrador[contador] - (3*strtoi(sequency[j+1]));
	j=j+3;
	}

valor=0;
contador=0;
j=1;
while (j <=length(sumatoria)) 
	{
	x=sumatoria[j]
	y=sumatoria[j+1]
	if(y==-2)
		{
		y=9;
		}
	if(y==-1)
		{
		y=8;
		}
	valor=strtoi(paste0(x,y))
	valor=valor+63;
	contador=contador+1;
	codificacion[contador]=chr(valor);
	j=j+2;
	}
valor=0;
valor2=0;
for (i in 1:length(sumatoria)) 
	{
	if(is.odd(i)==TRUE)
		{
		valor=valor+sumatoria[i];
		valor2=valor2+calibrador[i];
		}
	else
		{
		valor=valor-sumatoria[i];
		valor2=valor2-calibrador[i];
		}
	}

codificacion[19]=toString(valor);
codificacion[21]=toString(valor2);
return(codificacion);
}

#based on the first letter is the missing code (sum letter) ......
#for example if the code is 1000 we only will have 999, but the missin letter
#is on the first letter of the code.
write(seq2,"Datos.txt",100);
datos = scan("Datos_check.txt", what=character(1))
seq1=TwoWordSeq(datos);
splited = c(seq1[1],"-",seq1[2],"-",seq1[3],"-",seq1[4],"-",seq1[5],"-",seq1[6],"-",
seq1[7],"-",seq1[8],"-",seq1[9],"-",seq1[10],"-",seq1[11],"-",seq1[12],"-",
seq1[13],"-",seq1[14],"-",seq1[15],"-",seq1[16]);
write(c(datos[1],"-",splited),"compacted.txt",100,sep="");

contador = 0;
line=0;
ExtractoSeq = rep(0,102);
for (i in 1:length(datos)) 
	{
	contador=contador+1; 
	ExtractoSeq[contador]=datos[i];
	if ((contador==102)||(i==length(datos)))
		{
		line=line+1;
		contador=0;
		seq2=getNumSeq(ExtractoSeq);
		codificado=rep("0",17);
		codificado=codificador(seq2);
		write(c(codificado),"compacted.txt",1024,append = TRUE,sep="");
		ExtractoSeq = rep(0,102);
		}
	}

 



 
Decompressing Code in R

#library(ape);
#library(combinat);
library(stringr);
library(gtools)
#Program: DNA coding.
TransMatrix=rep(0,16);

is.even <- function(x) x %% 2 == 0
is.odd <- function(x) x %% 2 != 0

lectura = function(data)
{
	sub1=word(data[2],1,sep="!");
	sub2=word(data[2],2,sep="!");
	sub3=word(data[2],3,sep="!");
	decodificada=c(sub1,sub2,sub3);

for (i in 3:length(data)) 
	{
	sub1=word(data[i],1,sep="!")
	sub2=word(data[i],2,sep="!")
	sub3=word(data[i],3,sep="!")
	TransMatrix[i]=strtoi(substring);
	decodificada=c(decodificada,c(sub1,sub2,sub3))
	}
return(decodificada);
}
#######
numerifica = function(data)
{
largo=dim(data)[1];
ancho=dim(data)[2];
Transision=matrix(0,nrow=largo,ncol=36);
for (i in 1:largo)
	{
	for (j in 1:17)
		{
		valor=asc(data[i,j]);
		valor=strtoi(valor)-63;
		suma=sum(valor);
		if (length(valor)==1)
			{
			valor1=ifelse(valor>=0,valor%/%10,(valor*-1)%/%10); 
			valor2=ifelse(valor>=0,valor%%10,(valor*-1)%%10); 
			valor1=ifelse(valor>=0,valor1,(valor1*-1));	

			valor2=ifelse(valor2==8,-1,valor2);
			valor2=ifelse(valor2==9,-2,valor2);
			Transision[i,((j*2)-1)]=valor1;
			Transision[i,((j*2))]=valor2;
			}
		else
			{
			if(suma==197){valor=66}
			if(suma==209){valor=78}
			if(suma==211){valor=73}
			if(suma==217){valor=77}
			if(suma==218){valor=68}
			if(suma==231){valor=75}
			if(suma==319){valor=67}
			if(suma==323){valor=69}
			if(suma==325){valor=71}
			if(suma==326){valor=72}
			if(suma==331){valor=70}
			if(suma==339){valor=65}
			if(suma==341){valor=74}
			if(suma==350){valor=76}
			#print(c(probability[i,j],"suma:",sum(valor)))

			valor1=ifelse(valor>=0,valor%/%10,(valor*-1)%/%10); 
			valor2=ifelse(valor>=0,valor%%10,(valor*-1)%%10); 
			valor1=ifelse(valor>=0,valor1,(valor1*-1));	

			if(valor2==8){valor2=-1}
			if(valor2==9){valor2=-2}
			Transision[i,((j*2)-1)]=valor1;
			Transision[i,((j*2))]=valor2;
			}
		}

	Transision[i,35]=strtoi(data[i,18]);
	Transision[i,36]=strtoi(data[i,19]);	
	} 

return(Transision);
}
###
buscar = function(numero)
{
x=c(1,2,3,4);
valores =permutations(n=4,r=3,v=x,repeats.allowed=T);
tripleta2=rep(0,3);
suma=rep(0,64);
valores=cbind(valores,suma);
valores=cbind(valores,suma);
for (i in 1:64)
	{
	valores[i,4]=valores[i,1]+valores[i,2]-valores[i,3];
	valores[i,5]=((2*(valores[i,1]*valores[i,3]))-(3*valores[i,2]));
	} 

for (i in 1:64)
	{
	if(numero==valores[i,4])
		{
		tripleta1=c(valores[i,1],valores[i,2],valores[i,3]);
		###print(tripleta1)
		tripleta2=rbind(tripleta2,tripleta1);
		}
	}
largo=dim(tripleta2)[1];
x3=sample(2:largo,1);
cadenab=tripleta2[x3,];
return(cadenab);
}

####### transcript DNA

convierte_adn = function(data)
{
largo_ADN = sum(TransMatrix)+1;
cadena1=rep(0,34);
cadena2=rep(0,102);
cadena3=rep(0,34);
tripleta=rep(0,3);
secuencia=rep(0,1);
###Number of cicles for exact match

###real length chain "cadena" = dim(data)[1]
exacta=0;
distancia=100;
for (i in 1:1)
	{
	contador=0;
	cadena1 = data[i,1:34];
	check1=data[i,35]; 
	check2=data[i,36];
		for(repeticion in 1:cicles)
		{
		cadena2=buscar(cadena1[1]);
		cadena3[1]=((2*(cadena2[1]*cadena2[3]))-(3*cadena2[2]));

		for (j in 2:34)
			{		
			tripleta = buscar(cadena1[j]); 
			cadena2=c(cadena2,tripleta);
			cadena3[j]=((2*(tripleta[1]*tripleta[3]))-(3*tripleta[2]));
			}				
		###print(cadena1);
		###print(cadena2);
		###print(cadena3);
		suma1=0;
		suma2=0;
		for (k in 1:34)
			{
			suma1=ifelse(is.odd(k)==TRUE,suma1+cadena1[k],suma1-cadena1[k]); 
			suma2=ifelse(is.odd(k)==TRUE,suma2+cadena3[k],suma2-cadena3[k]);
			}
		exacta = ifelse(check2==suma2,1,0);
	
		if(exacta==1)
			{
			elegida=cadena2;
			repeticion=cicles;
			exacta==0;
			##print(elegida);
			##print(c(repeticion,suma1,suma2));
			##print(distancia);
			distancia=0;
			}
		else
			{
			if((abs(check2-suma2))<10)
				{
				elegida=cadena2;
				distancia=abs(check2-suma2);
				##print(elegida);
				##print(c(repeticion,suma1,suma2));
				##print(distancia);
				}
			}	
		
		}
	secuencia=c(secuencia,elegida);
	}

return(secuencia);
}
#######
recodifica= function(data)
{
text=substring(data[1], 1:17, 1:17)
text2=decodificada[c(2,3)];
DNA=c(text,text2);
i=4;
columnas=1;
while(i<=length(data))
	{
	text= substring(data[i], 1:17, 1:17)
	text2=decodificada[c(i+1,i+2)];
	text3=c(text,text2);
	DNA=rbind(DNA, text3);
	columnas=columnas+1;
	i=i+3;
	}
rownames(DNA)=c(1:columnas);
return(DNA);
}


codificador = function(sequency) 
{
contador=0;
calibrador=rep(0,34);
sumatoria=rep(0,34);
codificacion=rep("0",21);
codificacion[18]="/";
codificacion[20]="=";

j=1;
while (j<=length(sequency)) 
	{
	contador=contador+1;
	sumatoria[contador]=strtoi(sequency[j])+strtoi(sequency[j+1])-strtoi(sequency[j+2]);
	calibrador[contador]=strtoi(sequency[j])+strtoi(sequency[j+1])+strtoi(sequency[j+2]);
	j=j+3;
	}

valor=0;
contador=0;
j=1;
while (j <=length(sumatoria)) 
	{
	valor=strtoi(paste0(sumatoria[j],sumatoria[j+1]))
	valor=valor+33;
	contador=contador+1;
	codificacion[contador]=chr(valor);
	j=j+2;
	}
valor=0;
valor2=0;
for (i in 1:length(sumatoria)) 
	{
	if(is.odd(i)==TRUE)
		{
		valor=valor+sumatoria[i];
		valor2=valor2+calibrador[i];
		}
	else
		{
		valor=valor-sumatoria[i];
		valor2=valor2-calibrador[i];
		}
	}

codificacion[19]=toString(valor);
codificacion[21]=toString(valor2);
return(codificacion);
}


###measures the time machine for process the code
trials=1
tiempo <- proc.time()

#Main Decoding Sequence
datos = scan("compacted.txt", what=character(1));

for (i in 1:16) 
	{
	substring=word(datos[1],(i+1),sep="-")
	TransMatrix[i]=strtoi(substring);
	}
decodificada=lectura(datos);
probability = recodifica(decodificada);
numeros = numerifica(probability);
cicles=100000;
accuracy=80;
proc.time()-tiempo

for (x in 1:trials)
{
	t2 <- proc.time() ##function that measures time
	ADN_seq = convierte_adn(numeros);
	####ends time measuring sequence
	proc.time()-t2
	tiempo=rbind(tiempo,t2);
	accuracy=accuracy+10;
	cicles = x*100;
	tiempo[x+1,4]=cicles;
	tiempo[x+1,5]=accuracy;
}


