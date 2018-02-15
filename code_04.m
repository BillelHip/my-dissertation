load('matlab_02.mat');
delT = 0.1;

for loc = 1:1000

	%if(loc <= 3)

		MaxAgent = 100;
		Agents = [rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1)];
		for i = 2:MaxAgent
			Agents(i,:) = [rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1)];
		end

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

		LLa(1,loc) = 0.5;
		LLp(1,loc) = 0.5;
		LLv(1,loc) = 0.5;
		LLs(1,loc) = 0.5;

		for t = 1:numTime
			Ex(t) = 0.1;
			Er(t) = 0.1;
			Sc(t) = 0.1;
			Py(t) = 0.1;
			Tp(t) = 0.9;
			Abnorm(t) = 0.9;
			Jlbasic(t) = 0.9;
		end

			Ab(1) = aAb * Abnorm(1) + (1-aAb) * [La(1) * (1-Lp(1))] * Abnorm(1);
			Tc(1) = aTc * Ab(1) + (1-aTc) * Ex(1);
			Sr(1) = wSr1 * Py(1) + wSr2 * Er(1) + wSr3 * Sc(1); 

			Pr(1) = bPr * Sr(1) + (1-bPr) * Tc(1);
	
			Jl(1) = aJl * Jlbasic(1) + (1-aJl) * Lv(1);

			Jc(1) = [yJc * Ex(1) + (1-yJc) * Py(1)] * (1-Ls(1));
			Jd(1) = [uJd * Jl(1) + (1-uJd) * Tp(1)] * (1-Jc(1));
	

			Sv(1) = Tp(1) * (1-Jc(1));

	
			Pj(1) = Jd(1) * [1-(uPj * Jc(1) + (1-uPj) * Pr(1))];
			Aj(1) = [bAj* Pr(1)+(1-bAj)*Jc(1)]*(1-Jd(1));

			Js(1) = [yJs * Jd(1) + (1-yJs) * Ls(1)] * [1-[(wJs1 * Pr(1) + wJs2 * Er(1) + wJs3 * Jc(1)) * (1-Ls(1))]];

			Ss(1) = ySs * Js(1) + (1-ySs) * Sv(1);

			xTime = numTime/MaxAgent;
			LastxTime = xTime;
			PosAgent = 2;
	
			for t = 2:numTime
	
				whereTime = LastxTime+xTime;
				if(t >= whereTime)
					LastxTime = t;
					for i = 1:7
			
						AgentsMean = [Agents(1,i)];
						for j = 1:PosAgent
							AgentsMean(j,:) = [Agents(j,i)];
						end

						if(i == 1)	Ex(t) = mean(AgentsMean(:,1)); % max or mean
						elseif(i == 2)	Er(t) = mean(AgentsMean(:,1)); % max or mean
						elseif(i == 3)	Sc(t) = mean(AgentsMean(:,1)); % max or mean
						elseif(i == 4)	Py(t) = mean(AgentsMean(:,1)); % max or mean
						elseif(i == 5)	Tp(t) = mean(AgentsMean(:,1)); % max or mean
						elseif(i == 6)	Abnorm(t) = mean(AgentsMean(:,1)); % max or mean
						elseif(i == 7)	Jlbasic(t) = mean(AgentsMean(:,1)); % max or mean
						end
					end
					PosAgent = PosAgent + 1;
		
				end

				Ab(t) = aAb * Abnorm(t-1) + (1-aAb) * [La(t) * (1-Lp(t))] * Abnorm(t-1);
				Tc(t) = aTc * Ab(t) + (1-aTc) * Ex(t);
				Sr(t) = wSr1 * Py(t) + wSr2 * Er(t) + wSr3 * Sc(t); 

				Pr(t) = bPr * Sr(t) + (1-bPr) * Tc(t);
	
				Jl(t) = aJl * Jlbasic(t-1) + (1-aJl) * Lv(t);

				Jc(t) = [yJc * Ex(t) + (1-yJc) * Py(t)] * (1-Ls(t));
				Jd(t) = [uJd * Jl(t) + (1-uJd) * Tp(t)] * (1-Jc(t));
	

				Sv(t) = Tp(t) * (1-Jc(t));

	
				Pj(t) = Jd(t) * [1-(uPj * Jc(t) + (1-uPj) * Pr(t))];
				Aj(t) = [bAj* Pr(t)+(1-bAj)*Jc(t)]*(1-Jd(t));

				Js(t) = [yJs * Jd(t) + (1-yJs) * Ls(t)] * [1-[(wJs1 * Pr(t) + wJs2 * Er(t) + wJs3 * Jc(t)) * (1-Ls(t))]];

				Ss(t) = ySs * Js(t) + (1-ySs) * Sv(t);


				La(t) = La(t-1) + nLs * (Aj(t)-La(t-1)) * (1-La(t-1)) * La(t-1) * delT;
				Lp(t) = Lp(t-1) + nLp * (Pj(t)-Lp(t-1)) * (1-Lp(t-1)) * Lp(t-1) * delT;
				Lv(t) = Lv(t-1) + nLv * (Sv(t)-Lv(t-1)) * (1-Lv(t-1)) * Lv(t-1) * delT;
				Ls(t) = Ls(t-1) + nLs * (Ss(t)-Ls(t-1)) * (1-Ls(t-1)) * Ls(t-1) * delT;
		
				LLa(t,loc) = La(t);
				LLp(t,loc) = Lp(t);
				LLv(t,loc) = Lv(t);
				LLs(t,loc) = Ls(t);
	

			end % Time loop

	%else LLa(t,loc) = NaN; LLp(t,loc) = NaN; LLv(t,loc) = NaN; LLs(t,loc) = NaN;

	%end % if

end % Location loop

xlin=linspace(0,1000,1000);
ylin=linspace(0,1000,1000);
[X,Y]=meshgrid(xlin,ylin);
Z = griddata(xlin,ylin,LLa,X,Y,'cubic'); % LLp LLs LLv
mesh(X,Y,Z); 
