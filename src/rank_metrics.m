function [ errors ] = rank_metrics(calc_rank, range )
%rank_metrics Evaluate goodness of computed ranking

%metrics: correlation coefficient, sum of distances from true rank, sum of
%squared distances from true rank, max distance from true rank

abs_diff = sum(abs(calc_rank-(1:range)'));
square_diff = (calc_rank-(1:range)')'*(calc_rank-(1:range)');
max_diff = max(abs(calc_rank-(1:range)'));

%Get correlation coefficient
p= polyfit((1:range)',calc_rank,1);
yfit =  polyval(p,(1:range)');
yresid = calc_rank - yfit; 
SSresid = sum(yresid.^2);
SStotal = (length(calc_rank)-1) * var(calc_rank);
rsq = 1 - SSresid/SStotal;


errors = [abs_diff, square_diff, max_diff, rsq];
end

