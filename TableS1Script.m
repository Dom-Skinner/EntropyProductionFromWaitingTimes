%% Chlamydomonas data from Polin et al. Science, (2009).
A = readtable('Data/Chlamydomonas_science.txt');
t= A.Var1;
f_A = A.Var3/sum(A.Var3);

tA = sum(t.*f_A); % turning time
tA2 = sum(t.^2 .*f_A);
tB = 11.2; % run time from paper

sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For Chlamydomonas swim-turn dynamics, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /s\n',tA2/tA^2,sig)

%% Chlamydomonas data from Theves et al. Biophysical J., (2013).
A = readtable('Data/Theves2013.txt');
t= A.Var1;
f_A = A.Var2/sum(A.Var2);

tA = sum(t.*f_A); % run time
tA2 = sum(t.^2 .*f_A);
tB = 0.13; % turning time from paper
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For Pseudomonas putida swim-reversal time dynamics, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /s\n',tA2/tA^2,sig)

%% Myxobacteria data from Wu et al. PNAS, (2009).
A = readtable('Data/Wu2009.txt');
t= A.Var1;
f_A = A.Var2/sum(A.Var2);

tA = sum(t.*f_A)*60; % run time (in s)
tA2 = sum(t.^2 .*f_A) *3600;
tB = 0.0; % switching between states we assume is fast
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For Myxobacteria swarm-reversal time dynamics, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /s\n',tA2/tA^2,sig)

%% Paenibacillus dendritiformis data from Be'er et al. J. Bacteriology, (2013).
A = readtable('Data/Beer2013.txt');
t= A.Var1;
f_A = A.Var2/sum(A.Var2);

tA = sum(t.*f_A); % run time (in s)
tA2 = sum(t.^2 .*f_A);
tB = 4; % switching "a couple of seconds" from paper
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For Paenibacillus swarm-reversal time dynamics, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /s\n',tA2/tA^2,sig)

%% Myxococcus xanthus data from Sliusarenko et al. PNAS, (2006)
A = readtable('Data/Sliusarenko2006.txt');
t= A.Var1;
f_A = A.Var2/sum(A.Var2);

tA = sum(t.*f_A); % run time (in s)
tA2 = sum(t.^2 .*f_A);
tB = 0;
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For Myxococcus xanthus collective wave dynamics, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /s\n',tA2/tA^2,sig)

%% Cell repolarization time data from Zhou et al. PLOS one, (2020)
A = readtable('Data/Zhou2020.txt');
t= A.Var1;
f_A = A.Var2/sum(A.Var2);

tA = sum(t.*f_A)*60; % time in s
tA2 = sum(t.^2 .*f_A)*3600;
tB = 3*3600; % Total time period for oscillations ~10h ~ 2*4h + 2*tA
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For breast cancer cell repolarization, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.5f k_B /s\n',tA2/tA^2,sig)

%% Bi-stable perception switching from Krug et al. Proc. Royal Soc. B, (2008)
A = readtable('Data/Krug2008.txt');
t= A.Var1;
f_A = A.Var2/sum(A.Var2);

tA = sum(t.*f_A); % time in s
tA2 = sum(t.^2 .*f_A);
tB = tA; % A and B states equivalent
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For visual perception switching, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.5f k_B /s\n',tA2/tA^2,sig)

%% Seabird flight duration Elliott et al. Behavioral Ecology and Sociobiology, (2009)
A = readtable('Data/Elliott2009.txt');
t= A.Var1;
f_A = A.Var2/sum(A.Var2);

tA = sum(t.*f_A)*60; % time in s
tA2 = sum(t.^2 .*f_A)*3600;
tB = 8*60; % unclear exactly from paper, but around ~8mins between flights
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For Arctic seabird flight duration, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.5f k_B /s\n',tA2/tA^2,sig)

%% Bee flight duration Raine and Chittka Entomol. Gen. (2007)
A = readtable('Data/Raine2007.txt'); % (fig 5d LOT-LOT)
t= A.Var1;
f_A = A.Var2/sum(A.Var2);

tA = sum(t.*f_A); % time in s
tA2 = sum(t.^2 .*f_A);
tB = 5; % Fig 3, rough estimate
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For Bee flight duration, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.5f k_B /s\n',tA2/tA^2,sig)

%% Flagella rotation direction Bai et al. Science (2010)
tA = 20.5/1000; % read directly from fig 4(a) inset (middle).
tA2 = tA^2 + 357.8/1000^2; % they give the variance so convert
tB = 19.3/1000; 
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For Flagella switching, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.5f k_B /s\n',tA2/tA^2,sig)

%% Cell state switching Norman et al. Nature (2013)
A = readtable('Data/Norman2013.txt');
t= A.Var1;
f_A = A.Var2/sum(A.Var2);

tA = sum(t.*f_A)*60; % time in s
tA2 = sum(t.^2 .*f_A)*3600;

B = readtable('Data/Norman2013B.txt'); 
t= B.Var1;
f_B = B.Var2/sum(B.Var2);
tB = sum(t.*f_B)*60; % time in s

sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For cell state switching, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.5f k_B /s\n',tA2/tA^2,sig)

%% Binding-unbinding times for molecular motor Isojima et al. Nat. Chem. Biol. (2016)
A = readtable('Data/Isojima2016.txt');
t= A.Var1/100;
f_A = A.Var2/sum(A.Var2);

tA = sum(t.*f_A); % time in s
tA2 = sum(t.^2 .*f_A);

B = readtable('Data/Isojima2016B.txt'); % (fig 5d LOT-LOT)
t= B.Var1/100;
f_B = B.Var2/sum(B.Var2);
tB = sum(t.*f_B); % time in s

sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For molecular motor binding-unbinding, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.5f k_B /s\n',tA2/tA^2,sig)

%% Lifetime of RNA polymerase Cisse et al. Science (2013)
A = readtable('Data/Cisse2013.txt');
t= A.Var1;
f_A = A.Var2/sum(A.Var2);

tA = sum(t.*f_A); % time in s
tA2 = sum(t.^2 .*f_A);

tB = 500; % approximate time from paper

sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For RNA polymerase lifetime, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.5f k_B /s\n',tA2/tA^2,sig)