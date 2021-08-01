function tracadoparcial( XYZ )

%% Pack coordinates in list with streamlines separated by NaN.
   pl = reshape(permute(XYZ, [3,1,2]), [], 2);

   i = ~isnan(pl(:,1));
   j = i|[true;i(1:end-1)];
   pl = pl(j,:);

   % Pack streamline coordinates in a cell array suitable for use with
   % Matlab streamline, i.e., as in 'streamline(pollock(G, resSol));'
   flag = isnan(pl(:,1));
   ix = find(flag);
   dd  = diff([0;ix])-1;
   S = mat2cell(pl(~flag,:), dd, 2);
   
   Sl = reshape([S, repmat({[nan, nan]}, [numel(S),1])]',[], 1);

    Sl = vertcat(Sl{:});
 
    plot(Sl(:,1), Sl(:,2), 'r-');
   

end

