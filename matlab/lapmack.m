function [p,c]=lapmack(C)
% LAPMACK Linear Assignment Problem using Mack's approach -- Bradford 
% algorithm (as discussed by Jonker Volgenant)

% 
n=size(C,1);
dsum=0; % sum of slack values
dsol=zeros(1,n); % avoid modifying cost matrix
% step 1: Determine the basis in the cost matrix
[ci,bi]=min(C,[],2); % bi: base_i
% step 2a: Number of bases is each column 
cnt=accumarray(bi,1,[n,1]);
ba=zeros(n,1); % alternative basis, unpopulated
% step 3: select a column which contains more than one base
col=find(cnt>1,1); % first match suffices
% step 2b: Termination if every column contains one base
while(~isempty(col))
 row=find(bi==col);
 while true
%% step 4
 %% DBG
 % disp([bi,ci,cnt]')
 % D=C;  for i=1:n
 %     D(i,bi(i))=ci(i)+100;
 % %      E(i,bi(i))=E(i,bi(i))+100;
 % end 
 % D % DBG
 %% END
 %  for i=row'
 % set m = min { C_i,k-C_i,base_i | k (not?) in col}
 % determine delta=min (mi| i in ROW}
 % Let rr be a row and kk a column for which delta is assumed
  ncol=1:n;ncol(col)=[];
  [m,j]=min(C(row,ncol)+dsol(ncol),[],2);
  [d,r]=min(m-ci(row)); % d subtraction here
  rr=row(r);
  kk=ncol(j(r));

 % step 5: Audjust the dual solution
  dsol(col)=dsol(col)+d;
 % C(:,col)=C(:,col)+d;
  dsum=dsum+d*sum(col~=0);
  ci(row)=ci(row)+d;
  if cnt(kk)==0 % kk contains no base
    bi(rr)=kk;
    break; % out of step 4 loop
  else
 % step 7: column kk is base for some row
  % mark column kk as alternative base for row rr
   ba(rr)=kk;
  % col= col u {kk} row=row u { i| base_i=kk}
   col(end+1)=kk;
   bb=find(bi==kk);
   row(end+(1:length(bb)))=bb;
   %row(randperm(length(row)))=row; % perturb to avoid loop formation
  % goto step 4
  end
 end
  % Step 6: Augmentation

 % alter the current set of bases along the alternating path, starting in
 % column kk
% [bi,ba]'  %% DBG
 bi(rr)=kk; % enter directly
 bi(ba>0)=ba(ba>0);
 ba(:)=0;
 for i=1:n, j=bi(i);ci(i)=C(i,j)+dsol(j); end
 cnt=accumarray(bi,1,[n,1]);
 col=find(cnt>1,1); % prepare step 3
 % goto step 2
end
 %% DBG
 % D=C;  for i=1:n
 %     D(i,bi(i))=ci(i)+100;
 % end 
 % disp(D)  
 %% ENDDBG
%
c =-dsum; % cost: shift by sum of slack values
for i=1:n
    j=bi(i);
    c = c+C(i,j)+dsol(j); 
end
p=bi';
end

function test
%%
C=[3 7 6 6 ; 1 6 8 8 ; 3 0 8 1; 0 7 9 9]
[q,p]=lapmack(C)
[q,p]=lapmack(C')
[q,p]=lapmack([C,C;C,C])
[q,p]=lapmack([C,C;C,C]')

%%
end
