for i =1:NBETA
   hi = subplot(1,3,i);
   semilogy(alpha,squeeze(error_stats_zebrafish(2:5,3,:,i))')
   xlabel('\alpha');
   %ylabel('Rank correlation coefficient');
   %ylabel('Average absolute rank distance');
   %ylabel('Average rank distance squared');
   ylabel('Maximum rank distance');
   %title_str = sprintf('beta = %0.2f',beta(i)); 
   title_str = ['\beta =',num2str(beta(i))];
   title(title_str);
   legend('Pairwise','Fuzzy timestamps','Binary','Binary + affinities');
end
    