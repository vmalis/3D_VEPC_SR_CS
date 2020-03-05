function ANOVA=anovaCSStrain(a)

    fields=fieldnames(a);

    for i=1:size(fields,1)
 
        [p,~,stats] = anovan([a.(fields{i})],{[a.MVC] [a.Age]},'model',2,'varnames',{'MVC','Age'},'Display','off');
        
        posthoc=multcompare(stats,'Dimension',[1,2],'Display','off');
        ANOVA(1).(fields{i})= p(1);
        ANOVA(2).(fields{i})= p(2);
        ANOVA(3).(fields{i})= p(3);
        
        for j=1:size(posthoc,1)
            ANOVA(j+3).(fields{i})= posthoc(j,6);
         
        end
        
    end

end


