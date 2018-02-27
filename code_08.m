clc
clear
load('matlab_02.mat');
delT = 0.1;
numTime = 15000;

nLs = 0.3;
nLa = 0.3;
nLv = 0.3;
nLp = 0.3;

locNum = 5;
curAgent = 1;
for loc = 1:locNum

	%if(loc <= 3)

		MaxAgent = 10;
		%Agents = [rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1)];
        Agents = [((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1)];
        %[((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1)];
		for i = 1:MaxAgent
			%Agents(i,:) = [rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1)];

		Ex = zeros(1,numTime);
		Er = zeros(1,numTime);
		Sc = zeros(1,numTime);
		Py = zeros(1,numTime);
		Tp = zeros(1,numTime);
		Abnorm = zeros(1,numTime);
		Jlbasic = zeros(1,numTime);

		Ab = zeros(1,numTime);
		Tc = zeros(1,numTime);
		Pr = zeros(1,numTime);
		Sr = zeros(1,numTime);
		Jd = zeros(1,numTime);
		Jc = zeros(1,numTime);
		Jl = zeros(1,numTime);
		Sv = zeros(1,numTime);
		Pj = zeros(1,numTime);
		Aj = zeros(1,numTime);
		Js = zeros(1,numTime);
		Ss = zeros(1,numTime);
		La = zeros(1,numTime);
		Lp = zeros(1,numTime);
		Lv = zeros(1,numTime);
		Ls = zeros(1,numTime);

		La(1) = 0.5;
		Lp(1) = 0.5;
		Lv(1) = 0.5;
		Ls(1) = 0.5;
            for t = 1:numTime
                Ex(t) = Agents(i,1);
                Er(t) = Agents(i,2);
                Sc(t) = Agents(i,3);
                Py(t) = Agents(i,4);
                Tp(t) = Agents(i,5);
                Abnorm(t) = Agents(i,6);
                Jlbasic(t) = Agents(i,7);

                LEx(t,i) = Ex(t);
                LEr(t,i) = Er(t);
                LSc(t,i) = Sc(t);
                LPy(t,i) = Py(t);
                LTp(t,i) = Tp(t);
                LAbnorm(t,i) = Abnorm(t);
                LJlbasic(t,i) = Jlbasic(t);
            end
        
            Ab(1) = aAb * Abnorm(1) + (1-aAb) * [La(1) * (1-Lp(1))] * Abnorm(1);
            Tc(1) = aTc * Ab(1) + (1-aTc) * Ex(1);
            Sr(1) = wSr1 * Py(1) + wSr2 * Er(1) + wSr3 * Sc(1); 

            Pr(1) = bPr * Sr(1) + (1-bPr) * Tc(1);

            Jl(1) = aJl * Jlbasic(1) + (1-aJl) * Lv(1);

            Jc(1) = [yJc * Ex(1) + (1-yJc) * Py(1)] * (1-Ls(1));
            Jd(1) = [uJd * Jl(1) + (1-uJd) * Tp(1)] * (1-Jc(1));


            Sv(1) = Tp(1) * (1-Jc(1));


            Pj(1) = Pr(1) * [1-(uPj * Jc(1) + (1-uPj) * Jd(1))];
            Aj(1) = [bAj* Pr(1)+(1-bAj)*Jc(1)]*(1-Jd(1));

            Js(1) = [yJs * Jd(1) + (1-yJs) * Ls(1)] * [1-[(wJs1 * Pr(1) + wJs2 * Er(1) + wJs3 * Jc(1)) * (1-Ls(1))]];

            Ss(1) = ySs * Js(1) + (1-ySs) * Sv(1);
                
            %%%%%%%%%%%%%%%%%%%% Based on Time and Location 
            
            LAb(1,i) = Ab(1);
			LTc(1,i) = Tc(1);
			LSr(1,i) = Sr(1);

			LPr(1,i) = Pr(1);
	
			LJl(1,i) = Jl(1);

			LJc(1,i) = Jc(1);
			LJd(1,i) = Jd(1);

			LSv(1,i) = Sv(1);

			LPj(1,i) = Pj(1);
			LAj(1,i) = Aj(1);

			LJs(1,i) = Js(1);

			LSs(1,i) = Ss(1);
          
			for t = 2:numTime
                Ab(t) = aAb * Abnorm(t-1) + (1-aAb) * [La(t) * (1-Lp(t))] * Abnorm(t-1);
				Tc(t) = aTc * Ab(t) + (1-aTc) * Ex(t);
				Sr(t) = wSr1 * Py(t) + wSr2 * Er(t) + wSr3 * Sc(t); 

				Pr(t) = bPr * Sr(t) + (1-bPr) * Tc(t);
	
				Jl(t) = aJl * Jlbasic(t-1) + (1-aJl) * Lv(t);

				Jc(t) = [yJc * Ex(t) + (1-yJc) * Py(t)] * (1-Ls(t));
				Jd(t) = [uJd * Jl(t) + (1-uJd) * Tp(t)] * (1-Jc(t));

				Sv(t) = Tp(t) * (1-Jc(t));

				Pj(t) = Pr(t) * [1-(uPj * Jc(t) + (1-uPj) * Jd(t))];
				Aj(t) = [bAj* Pr(t)+(1-bAj)*Jc(t)]*(1-Jd(t));

				Js(t) = [yJs * Jd(t) + (1-yJs) * Ls(t)] * [1-[(wJs1 * Pr(t) + wJs2 * Er(t) + wJs3 * Jc(t)) * (1-Ls(t))]];

				Ss(t) = ySs * Js(t) + (1-ySs) * Sv(t);
                
                %%%%%%%%%%%%% Based on T , L

				LAb(t,i) = Ab(t);
				LTc(t,i) = Tc(t);
				LSr(t,i) = Sr(t);

				LPr(t,i) = Pr(t);
	
				LJl(t,i) = Jl(t);

				LJc(t,i) = Jc(t);
				LJd(t,i) = Jd(t);
	
				LSv(t,i) = Sv(t);

				LPj(t,i) = Pj(t);
				LAj(t,i) = Aj(t);

				LJs(t,i) = Js(t);

				LSs(t,i) = Ss(t);
                
				La(t) = La(t-1) + nLs * (Aj(t)-La(t-1)) * (1-La(t-1)) * La(t-1) * delT;
				Lp(t) = Lp(t-1) + nLp * (Pj(t)-Lp(t-1)) * (1-Lp(t-1)) * Lp(t-1) * delT;
				Lv(t) = Lv(t-1) + nLv * (Sv(t)-Lv(t-1)) * (1-Lv(t-1)) * Lv(t-1) * delT;
				Ls(t) = Ls(t-1) + nLs * (Ss(t)-Ls(t-1)) * (1-Ls(t-1)) * Ls(t-1) * delT;
		
				LLa(t,i) = La(t);
				LLp(t,i) = Lp(t);
				LLv(t,i) = Lv(t);
				LLs(t,i) = Ls(t);
                
                LLLa(t,i, loc) = LLa(t,i);
				LLLp(t,i, loc) = LLp(t,i);
				LLLv(t,i, loc) = LLv(t,i);
				LLLs(t,i, loc) = LLs(t,i);

			end % Time loop
            
            i = i + 1;
            Agents(i,:) = [((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1) ((0-1).*rand(1)+1)]; 
        end % Agent LOOP
        
	%else LLa(t,loc) = NaN; LLp(t,loc) = NaN; LLv(t,loc) = NaN; LLs(t,loc) = NaN;

	%end % if
    
    for t = 1:numTime
         meanLa(t,loc) = mean(LLa(numTime,:));
         meanLp(t,loc) = mean(LLp(numTime,:));
         meanLv(t,loc) = mean(LLv(numTime,:));
         meanLs(t,loc) = mean(LLs(numTime,:));
    end
    
end % Location loop
   

%%%%%%%%%%%%%%%%%%%%%% Stress long term
figure
xlin=linspace(0,numTime,numTime);
ylin=linspace(0,numTime,numTime);
[X,Y]=meshgrid(xlin,ylin);
%Z = griddata(xlin,ylin,LLs,X,Y,'natural'); % LLp LLs LLv
%v = X.^2 + Y.^3 - Z.^4;
axis tight;
%plot3(X,Y,LLs,'mo','MarkerSize',4);
%hold on
mesh(meanLs);
xlim([1 locNum]);
zlim([0 1]);
%mesh(X,Y,LLs);
%colorbar

%surfc(X,Y,LLs);

title('Graph of Stress long term');
xlabel('Location');
ylabel('Time');
zlabel('Level');
c = colorbar;
c.Label.String = 'Level of Effect';
set(c, 'ylim', [0 1]);

%%%%%%%%%%%%%%%%%%%%%% Overload long term
figure
xlin=linspace(0,numTime,numTime);
ylin=linspace(0,numTime,numTime);
[X,Y]=meshgrid(xlin,ylin);
%Z = griddata(xlin,ylin,LLv,X,Y,'natural'); % LLp LLs LLv

%axis tight;
%plot3(X,Y,LLv,'mo','MarkerSize',4);
%hold on

%mesh(X,Y,LLv);
mesh(meanLv);
%colorbar
axis tight;
xlim([1 locNum]);
zlim([0 1]);
%surfc(X,Y,LLv);

title('Graph of Overload long term');
xlabel('Location');
ylabel('Time');
zlabel('Level');
c = colorbar;
c.Label.String = 'Level of Effect';
set(c, 'ylim', [0 1]);

%%%%%%%%%%%%%%%%%%%%%% Passive Job long term
figure
xlin=linspace(0,numTime,numTime);
ylin=linspace(0,numTime,numTime);
[X,Y]=meshgrid(xlin,ylin);
%Z = griddata(xlin,ylin,LLp,X,Y,'natural'); % LLp LLs LLv

axis tight;
%plot3(X,Y,LLp,'mo','MarkerSize',4);
%hold on

%mesh(X,Y,LLp);
mesh(meanLp);
%colorbar
axis tight;
xlim([1 locNum]);
zlim([0 1]);
%surfc(X,Y,LLp);

title('Graph of Passive Job long term');
xlabel('Location');
ylabel('Time');
zlabel('Level');
c = colorbar;
c.Label.String = 'Level of Effect';
set(c, 'ylim', [0 1]);

%%%%%%%%%%%%%%%%%%%%%% Active Job long term
figure
xlin=linspace(0,numTime,numTime);
ylin=linspace(0,numTime,numTime);
[X,Y]=meshgrid(xlin,ylin);
%Z = griddata(xlin,ylin,LLa,X,Y,'natural'); % LLp LLs LLv


%plot3(X,Y,Z,'o','MarkerSize',4);
%hold on

%mesh(X,Y,LLa);

mesh(meanLa);
%colorbar;
axis tight;
xlim([1 locNum]);
zlim([0 1]);

title('Graph of Active Job long term');
xlabel('Location');
ylabel('Time');
zlabel('Level');
c = colorbar;
c.Label.String = 'Level of Effect';
set(c, 'ylim', [0 1]);