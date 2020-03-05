function ANOVA=anova1CSStrain(a)

    fields=fieldnames(a);

    for i=1:size(fields,1)

        [p,~,stats] = anova1([a.(fields{i})],[a.MVC],'off');

        posthoc=multcompare(stats,'Display','off');
        ANOVA(1).(fields{i})= p(1);
        %%ANOVA(2).(fields{i})= p(2);
        %%ANOVA(3).(fields{i})= p(3);
        size(posthoc)
        for j=1:size(posthoc,1)
            ANOVA(j+1).(fields{i})= posthoc(j,6);
        end
        
    end

end
