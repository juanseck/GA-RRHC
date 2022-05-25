%  Genetic Algorithm and Random-Restart Hill-Climbing (GA-RRHC)
%
%  Source codes demo version 1.0
%
%  Developed in MATLAB R2015a(7.08)
%
%  Author and programmer: Nayeli Jazmín Escamilla Serna, Juan Carlos Seck Tuoh Mora
%
% email:   es281355@uaeh.edu.mx
%          jseck@uaeh.edu.mx
%
%
%  Main paper:
%  Nayeli Jazmin Escamilla Serna, Juan Carlos Seck Tuoh Mora,
%  A hybrid search using genetic algorithms and hill-climbing with multiple restarts for Flexible Job Shop 
%  Scheduling Problems with high flexibility,
%  Journal of King Saud University - Computer and Information Sciences, 
%  DOI: http://
%_______________________________________________________________________________________________
% The initial parameters that you need are:
%__________________________________________
% numIndividuos = number of solutions
% numGeneraciones = maximum number of iterations
% numEstancamiento = maximum number of stagnation iterations
% numVecinos = number of neighbors for the CA-type nieghborhood to apply genetic operatos
% probElitista = probability of elitist solutions
% probMutacion = mutation probability
% iteracionesTotalesEscalada = number of iterations for the random restart hill-climbing (RRHC)
% iteracionesReinicioEscalada = number of iterations to restart the RRHC
% probOperCrit = probability of moving a critical operation in the RRHC
%______________________________________________________________________________________________

clear all
clc
rng shuffle

%Problem to solve
problema = 'mt20.fjs';

%GA-RRHC parameters
numIndividuos = 100;
numGeneraciones = 250;
numEstancamiento = 50;
numVecinos = 3;
probElitista = 0.02;
probMutacion = 0.1;
iteracionesTotalesEscalada = 100;
iteracionesReinicioEscalada = 30;
probOperCrit = 0.05;

%Read data problem
[numeroTrabajos, numeroMaquinas, numeroOperaciones, vectorNumOperaciones, vectorInicioOperaciones, vectorOperaciones, tablaTiempos, tablaMaquinasFactibles] = leerDatosProblema(problema);
%GA-RRHC algorithm
[mejorSO, mejorSM, mejorMakespan, PoblacionSO, PoblacionSM, PoblacionMakespan, convergencia, contIt] = GA_RRHC(numIndividuos, numGeneraciones, numEstancamiento, probElitista, numeroTrabajos, numeroMaquinas, numeroOperaciones, vectorNumOperaciones, vectorInicioOperaciones, vectorOperaciones, tablaTiempos, tablaMaquinasFactibles, numVecinos,probMutacion,iteracionesTotalesEscalada,iteracionesReinicioEscalada,probOperCrit,1);
%Best makespan calculated
disp(['Makespan: ' num2str(mejorMakespan)])
%Makespan convergence
figure(1)
plot(1:contIt,convergencia(1:contIt),'LineWidth',1.5)
ax = gca;
ax.FontSize = 13;
title('Makespan convergence')
xlabel('Iterations') 
ylabel('Makespan') 
xlim([1 contIt])
ylim([convergencia(end)-10 convergencia(1)+10])
grid on
%Gantt diagram
figure(2)
diagramaDeGanttMaquinas(mejorSO, mejorSM, numeroTrabajos, numeroMaquinas, numeroOperaciones, vectorInicioOperaciones, tablaTiempos)

%Read data
function [numeroTrabajos, numeroMaquinas, numeroOperaciones, vectorNumOperaciones, vectorInicioOperaciones, vectorOperaciones, tablaTiempos, tablaMaquinasFactibles] = leerDatosProblema(nombreArchivo)
%The function reads the data table of the problem and returns the number of jobs, the number of machines, 
%the vector with the number of operations per job, the vector with the previous operation number where each machine starts,
%the base vector with the number of operations per job, the base vector with the repeated operations for each job (1,...,1,2,...,2,...,n...,n...n),
%the table (operations/machines), the vector of the starting operation number of each job, the time table of each machine 
%and the feasible machines for each operation

%Data format:
%in the first line there are (at least) 2 numbers: the first is the number of jobs and the second the number
%of machines (the 3rd is not necessary, it is the average number of machines per operation)
%Every row represents one job: the first number is the number of operations of that job, the second number
%(let's say k>=1) is the number of machines that can process the first operation; then according to k, there are
%k pairs of numbers (machine,processing time) that specify which are the machines and the processing times;
%then the data for the second operation and so on...
%Example: Fisher and Thompson 6x6 instance, alternate name (mt06)
%6   6   1
%6   1   3   1   1   1   3   1   2   6   1   4   7   1   6   3   1   5   6
%6   1   2   8   1   3   5   1   5   10  1   6   10  1   1   10  1   4   4
%6   1   3   5   1   4   4   1   6   8   1   1   9   1   2   1   1   5   7
%6   1   2   5   1   1   5   1   3   5   1   4   3   1   5   8   1   6   9
%6   1   3   9   1   2   3   1   5   5   1   6   4   1   1   3   1   4   1
%6   1   2   3   1   4   3   1   6   9   1   1   10  1   5   4   1   3   1
%first row = 6 jobs and 6 machines 1 machine per operation
%second row: job 1 has 6 operations, the first operation can be processed by 1 machine that is machine 3 with processing time 1.
%Reference: http://people.idsia.ch/~monaldo/fjsp.html

%Open file
archivo=fopen(nombreArchivo,'r');
datos=fscanf(archivo,'%f');
numeroTrabajos=datos(1);
numeroMaquinas=datos(2);
vectorOperaciones=[];
%Indices to take jobs, operations and positions
indice=4;
nt=1;
%Loop for jobs
while(nt<=numeroTrabajos)
    vectorNumOperaciones(nt)=datos(indice);
    vectorInicioOperaciones(nt)=sum(vectorNumOperaciones(1:nt-1));
    operacionesTrabajo=ones(1,vectorNumOperaciones(nt))*nt;
    vectorOperaciones=[vectorOperaciones operacionesTrabajo];
    %Loop for operations
    for numOper=1:vectorNumOperaciones(nt)
        indice=indice+1;
        numMaq=datos(indice);
        %Loop for machines
        for i=1:numMaq
            indice=indice+1;
            maquina=datos(indice);
            indice=indice+1;
            tiempo=datos(indice);
            tablaTiempos(vectorInicioOperaciones(nt)+numOper,maquina)=tiempo;
        end
    end
    %Next job
    indice=indice+1;
    nt=nt+1;
end

%Operation number
numeroOperaciones=length(vectorOperaciones);
%Available machines per operation
tablaMaquinasFactibles=[];
for oper=1:length(tablaTiempos)
    indices_factibles = tablaTiempos(oper,:) ~= 0;
    tablaMaquinasFactibles=[tablaMaquinasFactibles; indices_factibles];
end
end

%Algorithm GA_RRHC
function [mejorSO, mejorSM, mejorMakespan, PoblacionSO, PoblacionSM, PoblacionMakespan, convergencia, contIt] = GA_RRHC(numIndividuos, numGeneraciones, numEstancamiento, probElitista, numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones, vectorInicioOperaciones, vectorOperaciones, tablaTiempos, tablaMaquinasFactibles,numVecinos,probMutacion,iteracionesTotalesEscalada,iteracionesReinicioEscalada,probOperCrit,bandImp)
%Initialize values
mejorSO = []; 
mejorSM = []; 
mejorMakespan = inf;
convergencia = [];
%Population vectors
PoblacionSO=zeros(numIndividuos,numOperaciones);
PoblacionSM=zeros(numIndividuos,numOperaciones);
PoblacionMakespan=zeros(numIndividuos,1);

%Table with the characteristics of each solution concerning the jobs,
%It is sorted by jobs and the order of their operations (J_11, J_12, ... Jnm-1, Jnm)
%Rows keep in this order the information:
%Machine assigned
%Processing position on the assigned machine
%End of operation time
%Operation duration
%Tail time 
%Operation position in SO
%Operation position in SM
PoblacionTablaTrabajos=zeros(6,numOperaciones,numIndividuos);

%Table with the characteristics of each solution concerning the machines,
%It is sorted by machines and the order of their operations (M_11, M_12, ... Mmo-1, Jmo)
%Rows keep in this order the information:
%Scheduled work
%Operation of scheduled work
%Final operation time
%Operation duration
%Tail time
%Operation position in SO
%Operation position in SM
PoblacionTablaMaquinas=zeros(6,numOperaciones,numIndividuos);

%Additional vectors to handle critical operations
PoblacionVectorMaquinas=zeros(numIndividuos,numeroMaquinas);
PoblacionVectorOrdenMaq=zeros(numIndividuos,numOperaciones);
PoblacionPosMk=zeros(numIndividuos,1);
PoblacionPosTT=zeros(numIndividuos,numOperaciones);
PoblacionPosTM=zeros(numIndividuos,numOperaciones);

%Loop to generate random population
for i=1:numIndividuos
    [PoblacionSO(i,:),PoblacionSM(i,:)] = generarIndividuoAleatorio(numeroMaquinas,numOperaciones,vectorOperaciones,tablaMaquinasFactibles);
end


%Evaluate population
[PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM] = calificarPoblacion(PoblacionSO, PoblacionSM, PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM, numeroTrabajos, numeroMaquinas, numOperaciones,numIndividuos,vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, 1);
%Best solution
[mejorSO, mejorSM, mejorMakespan] = mejorIndividuo(PoblacionSO, PoblacionSM, PoblacionMakespan);
%Number of elitist solutions
numIndEl=round(numIndividuos*probElitista);
if mod(numIndividuos-numIndEl,2)==1
    numIndEl = numIndEl+1;
end

%Iteration counter
contIt=1;
%Stagnation counter
contEst=1;
banderaCiclo=1;
%Convergence
convergencia=[];
convergencia(contIt)=mejorMakespan;

%Optimization loop
while(banderaCiclo)  
    %Selection
    [PoblacionSO, PoblacionSM, PoblacionMakespan,PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM] = seleccion(PoblacionSO, PoblacionSM, PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM, numIndividuos, numIndEl, 2);    
    %CA-like neihborhood
    [PoblacionSO,PoblacionSM,PoblacionMakespan,PoblacionPosMk,PoblacionTablaTrabajos,PoblacionTablaMaquinas,PoblacionVectorMaquinas,PoblacionVectorOrdenMaq,PoblacionPosTT,PoblacionPosTM] = vecindad(PoblacionSO, PoblacionSM, PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM, numIndividuos, numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, numVecinos, probMutacion, numIndEl);
    %RRHC
    [PoblacionSO,PoblacionSM,PoblacionMakespan,PoblacionPosMk,PoblacionTablaTrabajos,PoblacionTablaMaquinas,PoblacionVectorMaquinas,PoblacionVectorOrdenMaq,PoblacionPosTT,PoblacionPosTM] = busquedaPoblacionalEscaladaColina(PoblacionSO, PoblacionSM, PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM, numIndividuos, iteracionesTotalesEscalada, iteracionesReinicioEscalada, numeroTrabajos, numeroMaquinas, numOperaciones, vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones, numIndEl, probOperCrit);
    %Best new solution
    [nuevaMejorSO, nuevaMejorSM, nuevaMejorMakespan] = mejorIndividuo(PoblacionSO, PoblacionSM, PoblacionMakespan);
    %Update best solution if needed
    if nuevaMejorMakespan < mejorMakespan
        mejorSO = nuevaMejorSO;
        mejorSM = nuevaMejorSM;
        mejorMakespan = nuevaMejorMakespan;
        contEst=1;
    else
        %In other case, increase stagnation counter
        contEst=contEst+1;
    end
    %Halt condition
    if ((contEst>=numEstancamiento) || (contIt>=numGeneraciones))
        banderaCiclo=0;
    end
    %Print every 20 iterations
    if mod(contIt,20)==0 && bandImp==1
        disp(['Iteracion: ' num2str(contIt) ' Makespan: ' num2str(mejorMakespan)]) 
    end
    contIt=contIt+1;
    %Convergence record
    convergencia(contIt)=mejorMakespan; 
end
end

%Random smart-cell
function [so,sm] = generarIndividuoAleatorio(numeroMaquinas,numOperaciones,vectorOperaciones,tablaMaquinasFactibles)
%Operation permutation
indices=randperm(numOperaciones);
so=vectorOperaciones(indices);
%Random selection of feasible machines
lista_maquinas = 1:numeroMaquinas;
sm=zeros(1,numOperaciones);
for i=1:numOperaciones
    maquinas=tablaMaquinasFactibles(i,:);
    factibles = lista_maquinas(logical(maquinas));
    numMaquinas=length(factibles);
    indMaquina=randi([1 numMaquinas]);
    maquina=factibles(indMaquina);
    sm(i)=maquina;
end
end

%Makespan of all the population
function [PobMk, PobPMk, PobTT, PobTM, PobVM, PobVOM, PobPTT, PobPTM] = calificarPoblacion(PobSO, PobSM, PobMk, PobPMk, PobTT, PobTM, PobVM, PobVOM, PobPTT, PobPTM, numeroTrabajos, numeroMaquinas, numOperaciones,numIndividuos,vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, inicio)
for i=inicio:numIndividuos
    [PobMk(i), PobPMk(i), PobTT(:,:,i), PobTM(:,:,i), PobVM(i,:),PobVOM(i,:), PobPTT(i,:), PobPTM(i,:)] = calcularMakespanTablas(PobSO(i,:),PobSM(i,:),numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, PobTT(:,:,i), PobTM(:,:,i), PobVM(i,:), PobVOM(i,:), PobPTT(i,:), PobPTM(i,:));
end
end

%Makespan of each smart-cell
function [makespan,posmk,TablaTrabajos,TablaMaquinas,vectorMaquinas,vectorOrdenMaq,posTT,posTM] = calcularMakespanTablas(so,sm,numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM)
%Initialize makespan
makespan = 0;
posmk = 0;
%Final time of every machine
tiempoActualMaquina=zeros(1,numeroMaquinas);
%Current operation of every job
operActualTrabajo=ones(1,numeroTrabajos);
%Current time of every job
tiempoActualTrabajo=zeros(1,numeroTrabajos);
%Current operation of every machine
operActualMaquina=zeros(1,numeroMaquinas);

%Number of operations per machine in the smart-cell
for i=1:numeroMaquinas
    suma=sum(sm==i);
    vectorMaquinas(i)=suma;
end

%Operation loop
for i=1:numOperaciones
    trabajo=so(i);
    operacion=operActualTrabajo(trabajo);
    operActualTrabajo(trabajo)=operActualTrabajo(trabajo)+1;
    %Operation machine
    posicion=vectorInicioOperaciones(trabajo)+operacion;
    maquina=sm(posicion);
    %Update operation number in the machine
    operActualMaquina(maquina)=operActualMaquina(maquina)+1;
    operacionMaquina=operActualMaquina(maquina);
    %Processing time
    tiempo=tablaTiempos(posicion,maquina);
    %Initial time when the operation can start
    tiempoInicial=tiempoActualTrabajo(trabajo);
    aux=tiempoActualMaquina(maquina);
    if aux > tiempoInicial
        tiempoInicial=aux;
    end
    %Update times and makespan
    tiempoActualMaquina(maquina)=tiempoInicial+tiempo;
    tiempoActualTrabajo(trabajo)=tiempoInicial+tiempo;
    if makespan < tiempoInicial+tiempo
        makespan = tiempoInicial+tiempo;
        posmk = i;
    end
    %Update TablaTrabajos
    TablaTrabajos(1,posicion) = maquina; %machine assigned to the operation
    TablaTrabajos(2,posicion) = operacionMaquina; %position of the operation in the machine
    TablaTrabajos(3,posicion) = tiempo; %processing time of the operation
    TablaTrabajos(4,posicion) = tiempoInicial+tiempo; %final time of the operation
    TablaTrabajos(5,posicion) = 0; %tail time of the operation
    TablaTrabajos(6,posicion) = i; %operation position in SO
    %Update TablaMaquinas
    posicionMaquina=sum(vectorMaquinas(1:maquina-1))+operacionMaquina;
    TablaMaquinas(1,posicionMaquina) = trabajo; %job assigned to the machine
    TablaMaquinas(2,posicionMaquina) = operacion; %operation job
    TablaMaquinas(3,posicionMaquina) = tiempo; %processing time of the operation
    TablaMaquinas(4,posicionMaquina) = tiempoInicial+tiempo; %final time of the operation
    TablaMaquinas(5,posicionMaquina) = 0; %tail time of the operation
    TablaMaquinas(6,posicionMaquina) = i; %operation position in SO
    %Update other vectors
    vectorOrdenMaq(i) = maquina;
    posTT(i)=posicion;
    posTM(i)=posicionMaquina;
end

%Tail time calculus
for i=numOperaciones:-1:1
    postt=posTT(i);
    postm=posTM(i);
    maq=TablaTrabajos(1,postt);
    posm=TablaTrabajos(2,postt);
    trab=TablaMaquinas(1,postm);
    oper=TablaMaquinas(2,postm);
    if oper == vectorNumOperaciones(trab)
        tiempo_t=0;
    else
        tiempo_t=TablaTrabajos(5,postt+1)+TablaTrabajos(3,postt+1);
    end
    if posm==vectorMaquinas(maq)
        tiempo_m=0;
    else
        tiempo_m=TablaMaquinas(5,postm+1)+TablaMaquinas(3,postm+1);
    end
    tiempo_cola=max([tiempo_t tiempo_m]);
    TablaTrabajos(5,postt)=tiempo_cola;
    TablaMaquinas(5,postm)=tiempo_cola;
end
end

%Best smart-cell in the population
function [mejorSO, mejorSM, mejorMk, indice] = mejorIndividuo(PobSO, PobSM, PobMk)
[mejorMk, indice] = min(PobMk);
mejorSO=PobSO(indice,:);
mejorSM=PobSM(indice,:);
end

%Population selection by elitism and tournement 
function [nuevaPobSO, nuevaPobSM, nuevaPobMk, nuevaPobPosMk, nuevaPobTablaTrabajos, nuevaPobTablaMaquinas, nuevaPobVectorMaquinas, nuevaPobVectorOrdenMaq, nuevaPobPosTT, nuevaPobPosTM] = seleccion(PobSO, PobSM, PobMk, PobPosMk, PobTablaTrabajos, PobTablaMaquinas, PobVectorMaquinas, PobVectorOrdenMaq, PobPosTT, PobPosTM, numsoluciones, numElitista, numCompetidores)
%Elitism
[nuevaPobSO, nuevaPobSM, nuevaPobMk, nuevaPobPosMk, nuevaPobTablaTrabajos, nuevaPobTablaMaquinas, nuevaPobVectorMaquinas, nuevaPobVectorOrdenMaq, nuevaPobPosTT, nuevaPobPosTM] = seleccionElitista(PobSO, PobSM, PobMk, PobPosMk, PobTablaTrabajos, PobTablaMaquinas, PobVectorMaquinas, PobVectorOrdenMaq, PobPosTT, PobPosTM, numElitista);
ind=numElitista+1;
%Tournement
while (ind<=numsoluciones)
    indGanador = seleccionTorneo(PobMk, numsoluciones, numCompetidores);
    nuevaPobSO(ind,:)=PobSO(indGanador,:);
    nuevaPobSM(ind,:)=PobSM(indGanador,:);
    nuevaPobMk(ind)=PobMk(indGanador);
    nuevaPobPosMk(ind)=PobPosMk(indGanador);
    nuevaPobTablaTrabajos(:,:,ind)=PobTablaTrabajos(:,:,indGanador);
    nuevaPobTablaMaquinas(:,:,ind)=PobTablaMaquinas(:,:,indGanador);
    nuevaPobVectorMaquinas(ind,:)=PobVectorMaquinas(indGanador,:);
    nuevaPobVectorOrdenMaq(ind,:)=PobVectorOrdenMaq(indGanador,:);
    nuevaPobPosTT(ind,:)=PobPosTT(indGanador,:);
    nuevaPobPosTM(ind,:)=PobPosTM(indGanador,:);
    ind=ind+1;
end

end

%Elitist selection of numElitista individuals
function [nuevaPobSO, nuevaPobSM, nuevaPobMk, nuevaPobPosMk, nuevaPobTablaTrabajos, nuevaPobTablaMaquinas, nuevaPobVectorMaquinas, nuevaPobVectorOrdenMaq, nuevaPobPosTT, nuevaPobPosTM] = seleccionElitista(PobSO, PobSM, PobMk, PobPosMk, PobTablaTrabajos, PobTablaMaquinas, PobVectorMaquinas, PobVectorOrdenMaq, PobPosTT, PobPosTM, numElitista)
nuevaPobSO=PobSO;
nuevaPobSM=PobSM;
nuevaPobMk=PobMk;
nuevaPobPosMk=PobPosMk;
nuevaPobTablaTrabajos=PobTablaTrabajos;
nuevaPobTablaMaquinas=PobTablaMaquinas;
nuevaPobVectorMaquinas=PobVectorMaquinas;
nuevaPobVectorOrdenMaq=PobVectorOrdenMaq; 
nuevaPobPosTT=PobPosTT;
nuevaPobPosTM=PobPosTM;
[ ~, indices ] = sort( PobMk);
for i=1:numElitista
    nuevaPobSO(i,:)=PobSO(indices(i),:);
    nuevaPobSM(i,:)=PobSM(indices(i),:);
    nuevaPobMk(i)=PobMk(indices(i));
    nuevaPobPosMk(i)=PobPosMk(indices(i));
    nuevaPobTablaTrabajos(:,:,i)=PobTablaTrabajos(:,:,indices(i));
    nuevaPobTablaMaquinas(:,:,i)=PobTablaMaquinas(:,:,indices(i));
    nuevaPobVectorMaquinas(i,:)=PobVectorMaquinas(indices(i),:);
    nuevaPobVectorOrdenMaq(i,:)=PobVectorOrdenMaq(indices(i),:);
    nuevaPobPosTT(i,:)=PobPosTT(indices(i),:);
    nuevaPobPosTM(i,:)=PobPosTM(indices(i),:);
end
end

%Tournement selection among numCompetidores solutions
function indGanador= seleccionTorneo(PobMk,numsoluciones,numCompetidores)
indices=randperm(numsoluciones,numCompetidores);
[~,ind]=min(PobMk(indices));
indGanador=indices(ind);
end

%Apply genetic operators using a CA-like neighborhood
function [PoblacionSO, PoblacionSM, PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM] = vecindad(PoblacionSO, PoblacionSM, PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM, numIndividuos, numeroTrabajos, numeroMaquinas, numeroOperaciones, vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, numVecinos, probMutacion, numIndEl)
%Loop to take numVecinos for each solution
for i=1:numVecinos
    % OS Crossover taking 50% POX and 50% JBX, two-point crossover for MS
    % The first numElitista solutions are not considered
    % List of random crossovers
    listaC=randperm(numIndividuos-numIndEl)+numIndEl;
    listaC=reshape(listaC,floor((numIndividuos-numIndEl)/2),2);
    for j=1:floor((numIndividuos-numIndEl)/2)
        solucion1_so=PoblacionSO(listaC(j,1),:);
        solucion2_so=PoblacionSO(listaC(j,2),:);
        solucion1_sm=PoblacionSM(listaC(j,1),:);
        solucion2_sm=PoblacionSM(listaC(j,2),:);
        %OS crossover
        probCO=rand;
        if probCO<=0.5
            [solucion1_so,solucion2_so] = crucePOX(numeroTrabajos, solucion1_so, solucion2_so);
        else
            [solucion1_so,solucion2_so] = cruceJBX(numeroTrabajos, solucion1_so, solucion2_so);
        end
        [solucion1_sm,solucion2_sm] = cruceDosPuntos(numeroOperaciones,solucion1_sm,solucion2_sm);
        %SO Mutation, 50% swapping, 50% selecting 3 operations of different
        %jobs and swapping them
        probM=rand;
        if probM<=probMutacion
            probMO=rand;
            if probMO<=0.5
                solucion1_so = mutacionIntercambio(numeroOperaciones,solucion1_so);
            else
                solucion1_so = mutacionVecindad(numeroTrabajos,numeroOperaciones,solucion1_so);
            end
            solucion1_sm = mutacionMaquinas(numeroMaquinas,numeroOperaciones,solucion1_sm,tablaMaquinasFactibles);
        end
        probM=rand;
        if probM<=probMutacion
            probMO=rand;
            if probMO<=0.5
                solucion2_so = mutacionIntercambio(numeroOperaciones,solucion2_so);
            else
                solucion2_so = mutacionVecindad(numeroTrabajos,numeroOperaciones,solucion2_so);
            end
            solucion2_sm = mutacionMaquinas(numeroMaquinas,numeroOperaciones,solucion2_sm,tablaMaquinasFactibles);
        end
        %Evalute each solution and update if neccesary
        [mk, pmk, TT, TM, VM,VOM, PTT, PTM] = calcularMakespanTablas(solucion1_so,solucion1_sm,numeroTrabajos,numeroMaquinas,numeroOperaciones,vectorNumOperaciones,vectorInicioOperaciones,tablaTiempos,PoblacionTablaTrabajos(:,:,listaC(j,1)),PoblacionTablaMaquinas(:,:,listaC(j,1)),PoblacionVectorMaquinas(listaC(j,1),:),PoblacionVectorOrdenMaq(listaC(j,1),:),PoblacionPosTT(listaC(j,1),:),PoblacionPosTM(listaC(j,1),:));
        if mk<PoblacionMakespan(listaC(j,1))
            PoblacionSO(listaC(j,1),:)=solucion1_so;
            PoblacionSM(listaC(j,1),:)=solucion1_sm;
            PoblacionMakespan(listaC(j,1))=mk; 
            PoblacionPosMk(listaC(j,1))=pmk; 
            PoblacionTablaTrabajos(:,:,listaC(j,1))=TT;
            PoblacionTablaMaquinas(:,:,listaC(j,1))=TM;
            PoblacionVectorMaquinas(listaC(j,1),:)=VM;
            PoblacionVectorOrdenMaq(listaC(j,1),:)=VOM;
            PoblacionPosTT(listaC(j,1),:)=PTT;
            PoblacionPosTM(listaC(j,1),:)=PTM;
        end
        [mk, pmk, TT, TM, VM,VOM, PTT, PTM] = calcularMakespanTablas(solucion2_so,solucion2_sm,numeroTrabajos,numeroMaquinas,numeroOperaciones,vectorNumOperaciones,vectorInicioOperaciones,tablaTiempos,PoblacionTablaTrabajos(:,:,listaC(j,2)),PoblacionTablaMaquinas(:,:,listaC(j,2)),PoblacionVectorMaquinas(listaC(j,2),:),PoblacionVectorOrdenMaq(listaC(j,2),:),PoblacionPosTT(listaC(j,2),:),PoblacionPosTM(listaC(j,2),:));
        if mk<PoblacionMakespan(listaC(j,2))
            PoblacionSO(listaC(j,2),:)=solucion2_so;
            PoblacionSM(listaC(j,2),:)=solucion2_sm;
            PoblacionMakespan(listaC(j,2))=mk; 
            PoblacionPosMk(listaC(j,2))=pmk; 
            PoblacionTablaTrabajos(:,:,listaC(j,2))=TT;
            PoblacionTablaMaquinas(:,:,listaC(j,2))=TM;
            PoblacionVectorMaquinas(listaC(j,2),:)=VM;
            PoblacionVectorOrdenMaq(listaC(j,2),:)=VOM;
            PoblacionPosTT(listaC(j,2),:)=PTT;
            PoblacionPosTM(listaC(j,2),:)=PTM;
        end
    end
end
end

%Crossover POX
function [nuevaSol1,nuevaSol2] = crucePOX(numtrabajos, solucion1, solucion2)
%Crossover POX for sequences OS1 and OS2
%Split jobs in g1 and g2 at random
%Each element of OS1 is g1 is placed in OS1' at the same position
%the rest of places is defined by the elements in g2 at the same order in OS2
%OS2' is defined similarly, interchanging the action of g1 and g2
nuevaSol1= solucion1;
nuevaSol2= solucion2;
numg1=randi([1,numtrabajos-1]);
%Random permutation of jobs
per=randperm(numtrabajos);
g1=per(1:numg1);
g2=per(numg1+1:numtrabajos);
%OS1'
[~,indices1]=ismember(solucion1,g1);
indices1=indices1>0;
[~,indices2]=ismember(solucion2,g2);
indices2=indices2>0;
indices1N=(~indices1);
nuevaSol1(indices1N)=solucion2(indices2);
%OS2'
[~,indices1]=ismember(solucion2,g1);
indices1=indices1>0;
[~,indices2]=ismember(solucion1,g2);
indices2=indices2>0;
indices1N=(~indices1);
nuevaSol2(indices1N)=solucion1(indices2);
end

%Crossover JBX
function [nuevaSol1,nuevaSol2] = cruceJBX(numtrabajos, solucion1, solucion2)
%Crossover JBX for sequences OS1 and OS2
%Split jobs in g1 and g2 at random
%Each element of OS1 is g1 is placed in OS1' at the same position
%the rest of places is defined by the elements in g2 at the same order in OS2
%Each element of OS2 is g2 is placed in OS2' at the same position
%the rest of places is defined by the elements in g1 at the same order in
%OS1
nuevaSol1= solucion1;
nuevaSol2= solucion2;
numg1=randi([1,numtrabajos-1]);
%Random permutation of jobs
per=randperm(numtrabajos);
g1=per(1:numg1);
g2=per(numg1+1:numtrabajos);
[~,indices1]=ismember(solucion1,g1);
indices1=indices1>0;
[~,indices2]=ismember(solucion2,g2);
indices2=indices2>0;
indices1N=(~indices1);
nuevaSol1(indices1N)=solucion2(indices2);
%OS2'
indices2N=(~indices2);
nuevaSol2(indices2N)=solucion1(indices1);
end

%Two-point Crossover
function [nuevaSol1,nuevaSol2] = cruceDosPuntos(numOperaciones,solucion1,solucion2)
%For machine sequences
%Each sequence is divided in two random points
%Central sequences are exchanged to have new solutions
nuevaSol1= solucion1;
nuevaSol2= solucion2;
%Random cut places
per=randperm(numOperaciones);
lugar1=per(1);
if lugar1>per(2)
    lugar1=per(2);
    lugar2=per(1);
else
    lugar2=per(2);
end
%Crossover
nuevaSol1(1:lugar1)=solucion1(1:lugar1);
nuevaSol1(lugar1+1:lugar2)=solucion2(lugar1+1:lugar2);
nuevaSol1(lugar2+1:numOperaciones)=solucion1(lugar2+1:numOperaciones);
nuevaSol2(1:lugar1)=solucion2(1:lugar1);
nuevaSol2(lugar1+1:lugar2)=solucion1(lugar1+1:lugar2);
nuevaSol2(lugar2+1:numOperaciones)=solucion2(lugar2+1:numOperaciones);
end

%Swapping mutation
function [nuevaSol] = mutacionIntercambio(numOperaciones,solucion)
%Swap two random positions
nuevaSol = solucion;
pos=randperm(numOperaciones,2);
nuevaSol(pos(1))=solucion(pos(2));
nuevaSol(pos(2))=solucion(pos(1));
end

%Neighborhood mutation, swapping three operations from different jobs
function [nuevaSol] = mutacionVecindad(numTrabajos,numeroOperaciones,solucion)
nuevaSol = solucion;
posiciones=1:numeroOperaciones;
posicionesAl=zeros(1,3);
%Three random jobs
trab=randperm(numTrabajos,3);
%Random operations of each job
for i=1:3
    indices=solucion==trab(i);
    pos=posiciones(indices);
    posicionesAl(i)=pos(randi(length(pos),1));
end
%Permute positions
aux=randperm(3);
posicionesFin=posicionesAl(aux);
nuevaSol(posicionesFin)=solucion(posicionesAl);
end

%Machine mutation, change de feasible machine of half of the operations
function [nuevaSol] = mutacionMaquinas(numeroMaquinas, numeroOperaciones,solucion,tablaMaquinasFactibles)
nuevaSol = solucion;
lista_maquinas = 1:numeroMaquinas;
%Random selection of operations
indices=randperm(numeroOperaciones,floor(numeroOperaciones/2));
for i=1:length(indices)
    maquinas=tablaMaquinasFactibles(indices(i),:);
    factibles = lista_maquinas(logical(maquinas));
    numMaquinas=length(factibles);
    indMaquina=randi([1 numMaquinas]);
    maq=factibles(indMaquina);
    nuevaSol(indices(i))=maq;
end
end

%RRHC for all the population
function [PobSO, PobSM, PobMk, PobPMk, PobTT, PobTM, PobVM, PobVOM, PobPTT, PobPTM] = busquedaPoblacionalEscaladaColina(PobSO, PobSM, PobMk, PobPMk, PobTT, PobTM, PobVM, PobVOM, PobPTT, PobPTM, numIndividuos, iteracionesTotalesEscalada, iteracionesReinicioEscalada, numeroTrabajos, numeroMaquinas, numOperaciones,vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones, inicio, probOperCrit)
for i=inicio+1:numIndividuos
    [PobSO(i,:), PobSM(i,:), PobMk(i), PobPMk(i), PobTT(:,:,i), PobTM(:,:,i), PobVM(i,:), PobVOM(i,:), PobPTT(i,:), PobPTM(i,:)] = busquedaEscaladaColinaReinicioAleatorio(PobSO(i,:), PobSM(i,:), PobMk(i), PobPMk(i), PobTT(:,:,i), PobTM(:,:,i), PobVM(i,:), PobVOM(i,:), PobPTT(i,:), PobPTM(i,:), iteracionesTotalesEscalada, iteracionesReinicioEscalada, numeroTrabajos, numeroMaquinas, numOperaciones,vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones, probOperCrit);
end
end

%Random-restart hil-climbing over critical operations
function [mejorsolucionSO, mejorsolucionSM, mejorsolucionMk, mejorsolucionPMk, mejorsolucionTT, mejorsolucionTM, mejorsolucionVM, mejorsolucionVOM, mejorsolucionPTT, mejorsolucionPTM] = busquedaEscaladaColinaReinicioAleatorio(solucionSO, solucionSM, solucionMk, solucionPMk, solucionTT, solucionTM, solucionVM, solucionVOM, solucionPTT, solucionPTM, iteracionesTotalesEscalada, iteracionesReinicioEscalada, numeroTrabajos, numeroMaquinas, numOperaciones,vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones, probOperCrit)
mejorsolucionSO = solucionSO;
mejorsolucionSM = solucionSM;
mejorsolucionMk = solucionMk;
mejorsolucionPMk = solucionPMk;
mejorsolucionTT = solucionTT;
mejorsolucionTM = solucionTM;
mejorsolucionVM = solucionVM;
mejorsolucionVOM = solucionVOM;
mejorsolucionPTT = solucionPTT;
mejorsolucionPTM = solucionPTM;
%New solutions
nuevas_SO=zeros(iteracionesReinicioEscalada,numOperaciones);
nuevas_SM=zeros(iteracionesReinicioEscalada,numOperaciones);
%Loop for the RRHC
reinicio=0;
for i=1:iteracionesTotalesEscalada
    %In case of a new improved solution or a restart from the pile
    if reinicio==0
        %Critical operations and machines
        [RC, maquinasRC, tamRC] = rutaCriticaConTablas(numOperaciones, solucionTT, solucionTM, solucionPTT, solucionPTM, solucionMk, solucionPMk);
        indSMRC = RC(tamRC:-1:1);
        maqRC = maquinasRC(tamRC:-1:1);
        %Feasible machines for every critical operation
        maq_fac=zeros(tamRC,numeroMaquinas);
        num_maq_fac=zeros(tamRC,1);
        %Take the fesible machines different from teh current one
        lista_maquinas = 1:numeroMaquinas;
        for j=1:tamRC
            indMaqFac = tablaMaquinasFactibles(indSMRC(j),:);
            indMaqFac(maqRC(j))=false;
            factibles = lista_maquinas(logical(indMaqFac));
            num_maq_fac(j) = length(factibles);
            maq_fac(j,1:num_maq_fac(j))=factibles;
        end
        %Permutation of critical operations to avoid repetition in taking
        %the same operation
        indPermOC=randperm(tamRC);
        ind_inicial=0;
    end
    %Random selection of a critical operation and select a new feasible
    %machine
    smAux=solucionSM;
    ind_inicial=ind_inicial+1;
    if ind_inicial>tamRC
        ind_inicial=1;
    end
    pos_ruleta=indPermOC(ind_inicial);
    while num_maq_fac(pos_ruleta)==0
        ind_inicial=ind_inicial+1;
        if ind_inicial>tamRC
            ind_inicial=1;
        end
        pos_ruleta=indPermOC(ind_inicial);
    end
    pos_nueva=indSMRC(pos_ruleta);
    indice=randi([1, num_maq_fac(pos_ruleta)]);
    nueva_maquina=maq_fac(pos_ruleta,indice);
    smAux(pos_nueva)=nueva_maquina;
    soAux=solucionSO;
    bandera_oc=0;
    %Move a random critical operation with certain probability
    if rand<probOperCrit
        aux1=randi([1,tamRC]);
        aux2=indSMRC(aux1);
        pos_inicial=solucionTT(6,aux2);
        while(1)
            pos_final=randi([1,numOperaciones]);
            if pos_inicial ~= pos_final
                break
            end
        end
        %Swapping
        soAux(pos_inicial)=solucionSO(pos_final);
        soAux(pos_final)=solucionSO(pos_inicial);
        bandera_oc=1;
    end
    %Makespan of the new solution    
    if bandera_oc==1
        [makespan,posmk,TablaTrabajos,TablaMaquinas,vectorMaquinas,vectorOrdenMaq,posTT,posTM] = calcularMakespanTablas(soAux,smAux,numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, solucionTT, solucionTM, solucionVM, solucionVOM, solucionPTT, solucionPTM);
    else
        pos_os = solucionTT(6,pos_nueva);
        trabajo = soAux(pos_os);
        [makespan] = estimarMakespan(pos_os,trabajo,nueva_maquina,vectorNumOperaciones,tablaTiempos,solucionTT,solucionTM,solucionPTT,solucionPTM,solucionVM,solucionVOM);
    end
    %A beter solution has been found
    if makespan < mejorsolucionMk
        %Only if the makespan was estimated
        if bandera_oc==0
            [makespan,posmk,TablaTrabajos,TablaMaquinas,vectorMaquinas,vectorOrdenMaq,posTT,posTM] = calcularMakespanTablas(soAux,smAux,numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, solucionTT, solucionTM, solucionVM, solucionVOM, solucionPTT, solucionPTM);
        end
        %Update the best solution
        mejorsolucionSO = soAux;
        mejorsolucionSM = smAux;
        mejorsolucionMk = makespan;
        mejorsolucionPMk = posmk;
        mejorsolucionTT = TablaTrabajos;
        mejorsolucionTM = TablaMaquinas;
        mejorsolucionVM = vectorMaquinas;
        mejorsolucionVOM = vectorOrdenMaq;
        mejorsolucionPTT = posTT;
        mejorsolucionPTM = posTM;
        %Update the current solution
        solucionSO = soAux;
        solucionSM = smAux;
        solucionMk = makespan;
        solucionPMk = posmk;
        solucionTT = TablaTrabajos;
        solucionTM = TablaMaquinas;
        solucionVM = vectorMaquinas;
        solucionVOM = vectorOrdenMaq;
        solucionPTT = posTT;
        solucionPTM = posTM;
        reinicio=0; 
    else
        %It the estimated makespan is worst, increase the restart threshold
        reinicio = reinicio + 1;
        nuevas_SO(reinicio,:) = soAux;
        nuevas_SM(reinicio,:) = smAux;
    end
    %Restart de hill-climbing
    if reinicio >= iteracionesReinicioEscalada
        pos_al3=randi([1 reinicio]);
        soAux = nuevas_SO(pos_al3,:);
        smAux = nuevas_SM(pos_al3,:);
        [makespan,posmk,TablaTrabajos,TablaMaquinas,vectorMaquinas,vectorOrdenMaq,posTT,posTM] = calcularMakespanTablas(soAux,smAux,numeroTrabajos,numeroMaquinas,numOperaciones,vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, mejorsolucionTT, mejorsolucionTM, mejorsolucionVM, mejorsolucionVOM, mejorsolucionPTT, mejorsolucionPTM);
        solucionSO = soAux;
        solucionSM = smAux;
        solucionMk = makespan;
        solucionPMk = posmk;
        solucionTT = TablaTrabajos;
        solucionTM = TablaMaquinas;
        solucionVM = vectorMaquinas;
        solucionVOM = vectorOrdenMaq;
        solucionPTT = posTT;
        solucionPTM = posTM;
        reinicio = 0;
    end
end



end

%Critial path, going backward according to the previous operation
function [RC, maqRC, tamRC] = rutaCriticaConTablas(numOperaciones, TablaTT, TablaTM, PosTT, PosTM, Mk, PosMk)
RC=zeros(1,numOperaciones);
maqRC=zeros(1,numOperaciones);
tamRC=0;
tiempo_i = Mk;
ind = PosMk;
while(tiempo_i > 0)
    pos_t = PosTT(ind);
    pos_m = PosTM(ind);
    tiempo_i = tiempo_i - TablaTT(3,pos_t);
    tamRC = tamRC+1;
    RC(tamRC) = pos_t;
    maqRC(tamRC) = TablaTT(1,pos_t);
    if tiempo_i > 0
        if pos_t > 1
            if TablaTT(4,pos_t-1)==tiempo_i
                ind = TablaTT(6,pos_t-1);
            else
                ind = TablaTM(6,pos_m-1);
            end
        else
            ind = TablaTM(6,pos_m-1);
        end
    end
end
end

%Estimation of the makespan to speed computation 
function [nuevo_makespan] = estimarMakespan(pos_os,trabajo,nueva_maq,vectorNumOperaciones,tablaTiempos,solucionTT,solucionTM,solucionPTT,solucionPTM,solucionVM,solucionVOM)
    %Previous operations in the new machine
    ind_m=sum(solucionVOM(1:pos_os-1)==nueva_maq);
    ind_TM = sum(solucionVM(1:nueva_maq-1));
    %Final time of the previous operation in the new machine
    if ind_m == 0
        tf_ma = 0;
    else
        ind_TM = ind_TM + ind_m;
        tf_ma = solucionTM(4,ind_TM);
    end
    %Previous operations of the same job
    ind_o = solucionTM(2,solucionPTM(pos_os));
    ind_TT = solucionPTT(pos_os);
    %Final time of the previous operation in the same job
    if ind_o == 1
        tf_oa = 0;
    else
        tf_oa = solucionTT(4,ind_TT-1);
    end
    %Maximum previous final time
    max_tf_a = max([tf_ma, tf_oa]);
    %Time of the selected operation in the new machine
    duracion = tablaTiempos(ind_TT,nueva_maq);    
    %Procees time + Tile time of the following operation in the new machine
    if ind_m == solucionVM(nueva_maq)
        tc_ms = 0;
    else
        duracion_ms = solucionTM(3,ind_TM+1);
        tc_ms = duracion_ms + solucionTM(5,ind_TM+1);
    end
    %Procees time + Tile time of the following operation in the samejob
    if ind_o == vectorNumOperaciones(trabajo)
        tc_os = 0;
    else
        duracion_os = solucionTT(3,ind_TT+1);
        tc_os = duracion_os + solucionTT(5,ind_TT+1);
    end
    %Maximum following final time
    max_tf_s = max([tc_ms, tc_os]);
    %Estimated makespan
    nuevo_makespan = max_tf_a + duracion + max_tf_s;
end

%Draw the machine Gantt chart
function diagramaDeGanttMaquinas(so, sm, numeroTrabajos, numeroMaquinas, numOperaciones, vectorInicioOperaciones, tablaTiempos)
%First calculate makespan and keep a record of the operation and machine
%times
tiemposInMaq=zeros(numeroMaquinas,numeroTrabajos*2);
tiemposFinMaq=zeros(numeroMaquinas,numeroTrabajos*2);
trabMaq=zeros(numeroMaquinas,numeroTrabajos*2);
operTrabMaq=zeros(numeroMaquinas,numeroTrabajos*2);
makespan = 0;
tiempoActualMaquina=zeros(1,numeroMaquinas);
operActualMaquina=zeros(1,numeroMaquinas);
tiempoActualTrabajo=zeros(1,numeroTrabajos);
operActualTrabajo=zeros(1,numeroTrabajos);
for i=1:numOperaciones
    trabajo=so(i);
    operActualTrabajo(trabajo)=operActualTrabajo(trabajo)+1;
    operacion=operActualTrabajo(trabajo);
    posicion=vectorInicioOperaciones(trabajo)+operacion;
    maquina=sm(posicion);
    tiempo=tablaTiempos(posicion,maquina);
    tiempoInicial=tiempoActualTrabajo(trabajo);
    aux=tiempoActualMaquina(maquina);
    if aux > tiempoInicial
        tiempoInicial=aux;
    end
    tiempoActualMaquina(maquina)=tiempoInicial+tiempo;
    tiempoActualTrabajo(trabajo)=tiempoInicial+tiempo;
    operActualMaquina(maquina)=operActualMaquina(maquina)+1;
    if makespan < tiempoInicial+tiempo
        makespan = tiempoInicial+tiempo;
    end
    tiemposInMaq(maquina,operActualMaquina(maquina))=tiempoInicial;
    tiemposFinMaq(maquina,operActualMaquina(maquina))=tiempoInicial+tiempo;
    trabMaq(maquina,operActualMaquina(maquina))=trabajo;
    operTrabMaq(maquina,operActualMaquina(maquina))=operActualTrabajo(trabajo);
end

%Using the previous tables, the machine Gantt chart is displayed
colorT=colormap(hsv(numeroTrabajos));
colorM=colormap(lines(numeroMaquinas));
figure(2);
clf
grid on
grid minor
for maquina=1:numeroMaquinas
    y=[maquina-0.75 maquina-0.25 maquina-0.25 maquina-0.75];
    x=[-4 -4 -1 -1 ];
    patch(x,y,colorM(maquina,:));
    etiqueta=strcat( 'M', num2str(maquina));
    text(-3.5,maquina-0.5,etiqueta)
    y=[maquina-0.9 maquina-0.1 maquina-0.1 maquina-0.9];
    for oper=1:operActualMaquina(maquina)
        x = [tiemposInMaq(maquina,oper) tiemposInMaq(maquina,oper) tiemposFinMaq(maquina,oper) tiemposFinMaq(maquina,oper)];
        patch(x,y,colorT(trabMaq(maquina,oper),:));
        etiqueta=strcat( 'J', num2str(trabMaq(maquina,oper)), ',', num2str(operTrabMaq(maquina,oper)) ); 
        text(tiemposInMaq(maquina,oper),maquina-0.5,etiqueta)
    end
end
xlim([-4,makespan+5])
grafica=gca;
set(grafica,'YTick',[])
x=[0,0];
y=[0,numeroMaquinas];
hold on
plot(x,y,'k','LineWidth',1.25)
hold off
end









