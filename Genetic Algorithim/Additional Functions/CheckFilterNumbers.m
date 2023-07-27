[model_lin,model_loop1,model_loopN,model_brch] = model_generator_v5(5,'model_loopN');
totalModels=size(model_loopN,1)+size(model_loop1,1)+size(model_lin,1);
modelsWithCond=0;
state1A=2;
state2A=3;
state1B=1;
state2B=4;
for i=1:size(model_loopN,1)
    myBoolMat=(model_loopN{i,6})>0;
    if myBoolMat(state1A,state2A)>0 && myBoolMat(state1B,state2B)>0
        modelsWithCond=modelsWithCond+1;
    end
end
for i=1:size(model_loop1,1)
    myBoolMat=(model_loop1{i,6})>0;
    if myBoolMat(state1A,state2A)>0 && myBoolMat(state1B,state2B)>0
        modelsWithCond=modelsWithCond+1;
    end
end
for i=1:size(model_lin,1)
    myBoolMat=(model_lin{i,6})>0;
    if myBoolMat(state1A,state2A)>0 && myBoolMat(state1B,state2B)>0
        modelsWithCond=modelsWithCond+1;
    end
end