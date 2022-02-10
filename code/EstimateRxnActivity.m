load('..\..\data\CCLE\metabolic_data\Expression_Model_All');

load('C:/Projects/recon/Recon3D_301/Recon3D_301/Recon3DModel_301.mat'); 
model = Recon3DModel;
clear Recon3DModel

model.rules(find(ismember(model.rxns,'NADH2_u10mi'))) = {'(x(2220)&x(2216)&x(2198)&x(2207)&x(2202)&x(2224)&x(2201)&x(2190)&x(2217)&x(2189)&x(2206)&x(2183)&x(2204)&x(2218)&x(2211)&x(2195)&x(2222)&x(2203)&x(2219)&x(2185)&x(2186)&x(2213)&x(2227)&x(2215)&x(2210)&x(2180)&x(2197)&x(2194)&x(2188)&x(2182)&x(2181)&x(2200)&x(2205)&x(2214)&x(2196)&x(2187)&x(2192)&x(2223)&x(2226)&x(2209)&x(2228)&x(2184)&x(2212)&x(2191))'};
model.rules(find(ismember(model.rxns,'CYOOm3i'))) = {'((x(2232)|x(2247))&x(2241)&x(2231)&(x(2248)|x(2238))&(x(2246)|x(2242))&x(2239)&(x(2245)|x(2230)|x(2234))&(x(2237)|x(2244))&x(2235)&(x(2229)|x(2233))&x(2208)&x(2236)&x(2240)&x(2243))'};
model.rules(find(ismember(model.rxns,'ATPS4mi'))) = {'(x(2127)&x(2135)&x(2129)&x(2134)&x(2131))&((x(2137)&x(2138)&x(2139)&x(2152)&x(2156)&x(2158)&x(2143)&x(2144))|(x(2159)&x(2138)&x(2144)&x(2139)&x(2143)&x(2137)&x(2152)&x(2156))|(x(2138)&x(2144)&x(2139)&x(2143)&x(2137)&x(2157)&x(2152)&x(2156)))&(x(2160)&(x(2137)|x(2161))&x(2162))'};

% Estimating reaction expression for CCLE Cancer and Normal comparison
CT_ColNo = find(ismember(Expression_Model_All(1,:),'cancer_type'));
C_ColNo = find(ismember(Expression_Model_All(1,:),'cancer'));
S_ColNo = find(ismember(Expression_Model_All(1,:),'stroma'));
N_ColNo = find(ismember(Expression_Model_All(1,:),'normal'));
CCLE_ColNo = find(ismember(Expression_Model_All(1,:),'ccle_m'));
Gene_ColNo = find(ismember(Expression_Model_All(1,:),'ReconGeneName'));

cancer_type = unique(Expression_Model_All(2:end,CT_ColNo));
ExprReactions_All = [{'Reaction Abbreviation'},{'Subsystem'},{'CCLE'},{'cancer'},{'normal'},{'stroma'},{'cancer_type'}];

for i = 1:length(cancer_type)
    temp_exprDataCCLE = Expression_Model_All(find(ismember(Expression_Model_All(:,CT_ColNo),cancer_type(i))),[Gene_ColNo,CCLE_ColNo]);
    temp_exprDataCCLE = cat(1,temp_exprDataCCLE,cat(2,setdiff(model.genes,temp_exprDataCCLE(:,1)),num2cell(zeros(length(setdiff(model.genes,temp_exprDataCCLE(:,1))),1))));
    
    temp_exprDataCancer = Expression_Model_All(find(ismember(Expression_Model_All(:,CT_ColNo),cancer_type(i))),[Gene_ColNo,C_ColNo]);
    temp_exprDataCancer = cat(1,temp_exprDataCancer,cat(2,setdiff(model.genes,temp_exprDataCancer(:,1)),num2cell(zeros(length(setdiff(model.genes,temp_exprDataCancer(:,1))),1))));
    
    temp_exprDataNormal = Expression_Model_All(find(ismember(Expression_Model_All(:,CT_ColNo),cancer_type(i))),[Gene_ColNo,N_ColNo]);
    temp_exprDataNormal = cat(1,temp_exprDataNormal,cat(2,setdiff(model.genes,temp_exprDataNormal(:,1)),num2cell(zeros(length(setdiff(model.genes,temp_exprDataNormal(:,1))),1))));
    
    temp_exprDataStroma = Expression_Model_All(find(ismember(Expression_Model_All(:,CT_ColNo),cancer_type(i))),[Gene_ColNo,S_ColNo]);
    temp_exprDataStroma = cat(1,temp_exprDataStroma,cat(2,setdiff(model.genes,temp_exprDataStroma(:,1)),num2cell(zeros(length(setdiff(model.genes,temp_exprDataStroma(:,1))),1))));
    
    if (sum(cell2mat(temp_exprDataCCLE(:,2))) > 0)
        temp_expressionValCCLE = integrateExpr(model,temp_exprDataCCLE);
    else
        temp_expressionValCCLE = repmat(0,length(model.rxns),1);
    end
    temp_expressionValCancer = integrateExpr(model,temp_exprDataCancer);
    temp_expressionValStroma = integrateExpr(model,temp_exprDataStroma);
    temp_expressionValNormal = integrateExpr(model,temp_exprDataNormal);
    ExprReactions_All = cat(1,ExprReactions_All,cat(2,model.rxns,model.subSystems, num2cell(temp_expressionValCCLE),num2cell(temp_expressionValCancer),...
        num2cell(temp_expressionValNormal), num2cell(temp_expressionValStroma),...
        repmat(cancer_type(i),length(model.rxns),1)));
    disp(i)
end
ExprReactions_All_nonzero = [ExprReactions_All(1,:);ExprReactions_All(find(sum(cell2mat(ExprReactions_All(2:end,[3:6])),2))+1,:)];

% save results
% save('..\..\data\CCLE\metabolic_data\min_max_new_model\ExprReactions_All_nonzero','ExprReactions_All_nonzero')
% save('..\..\data\CCLE\metabolic_data\min_max_new_model\ExprReactions_All','ExprReactions_All')
% writetable(cell2table(ExprReactions_All_nonzero),'..\..\data\CCLE\metabolic_data\min_max_new_model\ExprReactions_All_nonzero.dat','Delimiter',',','WriteVariableNames',false)
% writetable(cell2table(ExprReactions_All),'..\..\data\CCLE\metabolic_data\min_max_new_model\ExprReactions_All.dat','Delimiter',',','WriteVariableNames',false)
