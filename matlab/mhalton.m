function y=mhalton(n,k)
% MHALTON Mingled Halton Sequence of Manal Bayousef

if(k>100), error('Out of multipliers'); end
y=zeros(n,k);
% Data from the thesis:
% prime numbers and optimal multipliers for dimensions up to 100
vec=@(x)x(:)';
PM=[ ...
2 1 73 42 179 50 283 104 419 302
3 2 79 60 181 140 293 105 421 88
5 3 83 24 191 74 307 173 431 93
7 5 89 54 193 149 311 115 433 307
11 7 97 56 197 52 313 122 439 136
13 11 101 86 199 129 317 84 443 312
17 14 103 84 211 39 331 136 449 103
19 15 107 59 223 48 337 70 457 140
23 14 109 85 227 178 347 92 461 101
29 21 113 24 229 188 349 62 463 176
31 13 127 97 233 84 353 149 467 181
37 22 131 76 239 70 359 126 479 145
41 15 137 52 241 55 367 97 487 176
43 12 139 98 251 98 373 99 491 222
47 10 149 115 257 69 379 274 499 135
53 39 151 111 263 71 383 157 503 195
59 34 157 114 269 194 389 105 509 213
61 43 163 101 271 209 397 155 521 380
67 18 167 35 277 78 401 327 523 199
71 21 173 76 281 177 409 97 541 152];
P=vec(PM(:,1:2:end)); % prime
M=vec(PM(:,2:2:end)); % multiplier
R=P(1:k);
S=M(1:k);
for i=1:n
    j=i*ones(1,k);
    Q=1./R;
    while(any(j))
        l=mod(j,R);
        y(i,:)=y(i,:)+mod(S.*l,R).*Q;
        j=(j-l)./R;
        Q=Q./R;
    end
end
end

function y=myhalton(n,k)
%MYHALTON Low Discrepancy Series
y=zeros(n,k);
% primes
P=[2,3,5,7,11,13,17,19,23,29,31,37]; % dim>=k
R=P(1:k);
for i=1:n
    j=i*ones(1,k);
    Q=1./R;
    while(any(j))
        l=mod(j,R);
        y(i,:)=y(i,:)+l.*Q;
        j=(j-l)./R;
        Q=Q./R;
    end
end
end

function testmyhalton
%%
M=12;
v=myhalton(2000,M);
z=mhalton(2000,M);
for i=1:M-1;for j=(i+1):M
        subplot(1,2,1)
        plot(z(:,i),z(:,j),'.');
          subplot(1,2,2)
        plot(v(:,i),v(:,j),'.');
        pause(2);
end;end
%%
end