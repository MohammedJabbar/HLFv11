from fractions import Fraction as frac
from operator import add
from operator import neg
from operator import mod
import re
import time
import binascii
from threading import Thread
import Compression_Ratio


def Bin2Poly(message):
  global N
  res2 = ''.join([decimalToBinary(ord(i)) for i in message])
  print(len(res2))
  #message = [ord(c) for c in message]
  #message = [ str('0'*(10-len(decimalToBinary(D))))+ decimalToBinary(D)+'1' for D in message]
  time_start = time.process_time()
  message = [res2[i-N:i] for i in range(N,len(res2)+1,N)]
  print('message len :',len(message))
  # message = [ Encrypted_mas(i) for i in message ]
  message = [ Compression_Ratio.encoderLen(Encrypted_mas(i)) for i in message ]
  #print(len(res2),len(message[0][0]))
  return str(message),(time.process_time()- time_start)
  
def Bin2Poly2(message):
  global N
  res2 = ''.join([decimalToBinary(ord(i)) for i in message])
  #message = [ord(c) for c in message]
  #message = [ str('0'*(10-len(decimalToBinary(D))))+ decimalToBinary(D)+'1' for D in message]
  time_start = time.process_time()
  try:
    message = [res2[i-N:i+1] for i in range(N,len(res2)+1,N)]
  except:
    message.append(res2[(-1*N):])
  # message = [ Encrypted_mas(i) for i in message ]
  message = [ Encrypted_mas(i) for i in message ]
  #print(len(res2),len(message[0][0]))
  return str(message),(time.process_time()- time_start)


# def TrimBin(message):
#    NewMessage = []
#    for n,i in enumerate(message):
#      cont = 0
#      for j in i:
#        if j == '1':
#          break
#        cont += 1
#      NewMessage.append(message[n][cont:][::-1])
#    return NewMessage

def modPoly(c,k):
  if(k==0):
    pass
    # print("Error in modPoly(c,k). Integer k must be non-zero")
  else:
    return [fracMod(x,k) for x in c]

def subPoly(c1,c2):
  [c1,c2]=resize(c1,c2)	
  c2=list(map(neg,c2))
  out=list(map(add, c1, c2))
  return trim(out)

def multPoly(c1,c2):
  order=(len(c1)-1+len(c2)-1)
  out=[0]*(order+1)
  for i in range(0,len(c1)):
    for j in range(0,len(c2)):
      out[j+i]=out[j+i]+c1[i]*c2[j]
  return trim(out)

def resize(c1,c2):
  if(len(c1)>len(c2)):
    c2=c2+[0]*(len(c1)-len(c2))
  if(len(c1)<len(c2)):
    c1=c1+[0]*(len(c2)-len(c1))
  return [c1,c2]

def trim(seq):
  if len(seq) == 0:
    return seq
  else:
    for i in range(len(seq) - 1, -1, -1):
      if seq[i] != 0:
        break
    return seq[0:i+1]

def extEuclidPoly(a,b):
  switch = False
  a=trim(a)
  b=trim(b)
  if len(a)>=len(b):
    a1, b1 = a, b
  else:
    a1, b1 = b, a
    switch = True
  Q,R=[],[]
  while b1 != [0]:
    [q,r]=divPoly(a1,b1)
    Q.append(q)
    R.append(r)
    a1=b1
    b1=r
  S=[0]*(len(Q)+2)
  T=[0]*(len(Q)+2)

  S[0],S[1],T[0],T[1] = [1],[0],[0],[1]

  for x in range(2, len(S)):
    S[x]=subPoly(S[x-2],multPoly(Q[x-2],S[x-1]))
    T[x]=subPoly(T[x-2],multPoly(Q[x-2],T[x-1]))

  gcdVal=R[len(R)-2]
  s_out=S[len(S)-2]
  t_out=T[len(T)-2]
	### ADDITIONAL STEPS TO SCALE GCD SUCH THAT LEADING TERM AS COEF OF 1:
  scaleFactor=gcdVal[len(gcdVal)-1]
  gcdVal=[x/scaleFactor for x in gcdVal]
  s_out=[x/scaleFactor for x in s_out]
  t_out=[x/scaleFactor for x in t_out]
	
  if switch:
    return [gcdVal,t_out,s_out]
  else:
    return [gcdVal,s_out,t_out]

def divPoly(N,D):
	N, D = list(map(frac,trim(N))), list(map(frac,trim(D)))
	degN, degD = len(N)-1, len(D)-1
	if(degN>=degD):
		q=[0]*(degN-degD+1)
		while(degN>=degD and N!=[0]):
			d=list(D)
			[d.insert(0,frac(0,1)) for i in range(degN-degD)]
			q[degN-degD]=N[degN]/d[len(d)-1]
			d=[x*q[degN-degD] for x in d]
			N=subPoly(N,d)
			degN=len(N)-1
		r=N	
	else:
		q=[0]
		r=N
	return [trim(q),trim(r)]

def addPoly(c1,c2):
	[c1,c2]=resize(c1,c2)
	out=list(map(add, c1, c2))
	return trim(out)

def cenPoly(c,q):
	u=float(q)/float(2)
	l=-u
	c=modPoly(c,q)
	c=[mod(x,-q) if x>u else x for x in c]
	c=[mod(x,q) if x<=l else x for x in c]
	return c

def reModulo(num,div,modby):
	[_,remain]=divPoly(num,div)
	return modPoly(remain,modby)

def cn(a, b):
    x,y, u,v = 0,1, 1,0
    while a != 0:
        q, r = b//a, b % a
        m, n = x-u * q, y-v * q
        b,a, x,y, u,v = a,r, u,v, m,n
    return b ,x,y

def egcd(a, b):
    b,x,y = cn(a, b)
    gcdVal = b
    return gcdVal, x, y


#Modular inverse
#An application of extended GCD algorithm to finding modular inverses:
def modinv(a, m):
    gcdVal, x, y = egcd(a, m)
    if gcdVal != 1:
        return None  # modular inverse does not exist
    else:
        return x % m

#Modulus Function which handles Fractions aswell
def fracMod(f,m):
  [tmp,t1,t2]=egcd(f.denominator,m)
  if tmp!=1:
    # print ("ERROR GCD of denominator and m is not 1")
    1/0
  else:
    out=modinv(f.denominator,m)*f.numerator % m
  return out

def NT(Np=167,pp=3,qp=128,fp=[-1,1,1,0,-1,0,1,0,0,1,-1],gp=[-1,0,1,1,0,1,0,0,-1,0,-1]):
  global h,D,f_p,f_q,randPol,N,p,q,f,g,text,ciphertext
  N=Np
  p=pp
  q=qp
  f=fp
  g=gp 
  text="None"
  ciphertext ="None"

  # print("==== Bob generates public key =====")
  # print("Values used:")
  # print(" N=",N)
  # print(" p=",p)
  # print(" q=",q)
  # print("========")
  # print("\nBob picks two polynomials (g and f):")

  f=[-1,0,1,1,-1,0,-1]
  g=[0,-1,-1,0,1,0,1]

  d=2

  # print("f(x)= ",f)
  # print("g(x)= ",g)

  D=[0]*(N+1)
  D[0]=-1
  D[N]=1

  # print("\n====Now we determine F_p and F_q ===")
  [gcd_f,s_f,t_f]=extEuclidPoly(f,D)

  f_p=modPoly(s_f,p)
  f_q=modPoly(s_f,q)
  # print("F_p:",f_p)
  # print("F_q:",f_q)

  x=multPoly(f_q,g)
  h=reModulo(x,D,q)

  # print("\n====And finally h====")
  # print("f_q x g: ",x)
  # print("H (Bob's Public Key): ",h)

  # print("\n====Let's encrypt====")
  randPol=[-1,-1,1,1]
  return [h,N,p,q]

# print("Random:\t\t\t",randPol)
def Encrypted_mas(msg):
    msg = [int(ii) for ii in msg]
    e_tilda=addPoly(multPoly(multPoly([p],randPol),h),msg)
    time_start1 = time.process_time()
    e=reModulo(e_tilda,D,q)
    return e

# print("\n====Let's decrypt====")

def decrypt_mas(e):
    time_start1 = time.process_time()
    time_start2 = time.process_time()
    tmp=reModulo(multPoly(f,e),D,q)
    print("time-1:",time.process_time()- time_start1)
    
    time_start1 = time.process_time()
    centered=cenPoly(tmp,q)
    print("time-2:",time.process_time()- time_start1)
    
    time_start1 = time.process_time()
    m1=multPoly(f_p,centered)
    print("time-3:",time.process_time()- time_start1)
    
    time_start1 = time.process_time()
    tmp=reModulo(m1,D,p)
    print("time-4:",time.process_time()- time_start1)
    print("time-all:",time.process_time()- time_start2)
  
    return tmp

def decimalToBinary(n):  
    return bin(n).replace("0b", "")

def splitData(x):
    junkers = re.compile('[[" \]]')
    result = junkers.sub('', x).split(',')
    # print('reslt',(result))
    s = []
    s2 = []
    #print('len reslt ',len(result[0]))
    for i in result :
      try:
        s2.append(int(i))
      except:
         pass
        # print(x)
      if len(s2) == N:
          s.append(s2.copy())
          s2=[]
    if len(s2) < N and len(s2):
        s.append(s2.copy())
    return s
    
def listData2Msg(s):
  #print("s len. :",len(s))
  #print(s)
  message = [ ''.join(str(e) for e in trim(decrypt_mas(s[i]))) for i in range(0,len(s)) ]
  #print(len(message))
  #print(len(message[0]))
  listChar = ''
  for j in range(len(message)):
    for i in  range(7,len(message[j])+1,7):
      try:
        listChar += chr(int(''.join(message[0][i-7:i]),2))
      except:
        listChar += chr(int(''.join(message[0][(i-7)*-1:]),2))
  return listChar

def DEC(s):
  f = [decrypt_mas(i) for i in s]
  return f

def TestCode():
  for lk in zip([251 ,401 ,439,487,593,743],[2048,2048,2048,2048,2048,2048],[2008,3033,3501,4383,5193,7690]):
    pa = NT(lk[0],3,lk[1])
    time_end = []
    for A in [lk[2]]:
      byte = int(A/8)
      byte = int(byte/3)
      my_str = "Z"*byte
      #my_str = "1234567890123456"
      lenS = len(my_str)*8
      #my_str_as_bytes = str.encode(my_str)
      #my_str_as_bytes = b'&\x1cu\x8e\xa2\x97\x195\xba\x86<#\xfc9\x9b\x1c'
      #print('len s:',len(my_str))
      x ,t= Bin2Poly(my_str)
      s = splitData(x)
      s = [Compression_Ratio.decoderLen(i) for i in (s)]
      #print("s len :",len(s))
      #print("s Type :",type(s))
      #s = [decrypt_mas(i) for i in s]
      #listChar = listData2Msg(s)
      #print(f">> {listChar}")
      n = 1
      for N in range(1,n+1):
        time_start = time.process_time()
        #threads = []
        #for _ in range(N) : # range(multiprocessing.cpu_count())
          #threads.append(Thread(target = Bin2Poly, args = (my_str,)))
        #for thread in threads: 
          #thread.start()
        #for thread in threads: 
          #thread.join()
          # time_end.append(time.process_time()- time_start)
        print(len(s))
        f = [decrypt_mas(i) for i in s]
        Time = time.process_time()- time_start
        print(f'Time : {Time} ,cycles : {Time*1400_000_000} , of {N} message , bit-{lenS} , N:{pa[1]} , q:{pa[-1]}')

if __name__ == '__main__':
   TestCode()
  # for lk in zip([167,251,347,503],[128,128,256,256]):
  #   pa = NT(lk[0],3,lk[1])
    # time_end = []
    # for A in [16,24,32]:
    #   byte = A
    #   my_str = "Z"*byte
    #   lenS = len(my_str)*8
    #   # my_str_as_bytes = str.encode(my_str)
    #   x = str(Bin2Poly(my_str))
    #   s = splitData(x)
    #   s = [decrypt_mas(i) for i in s]
      
    #   listChar = listData2Msg(s)
    #   print(f">> {listChar}")
    #   s = splitData(x)
    #   n = 25
    #   for N in range(n):
    #     time_end = []
    #     for _ in range(N+1):
    #           time_start = time.process_time()
    #           s = [decrypt_mas(i) for i in s]
    #           time_end.append(time.process_time()- time_start)
    #           listChar = listData2Msg(s)
    #     Time = sum(time_end)/1
    #     print(f'Time : {Time} , of {N+1} message , AES-{lenS} , N:{pa[1]} , q:{pa[-1]}')
