exc=0;
act=zeros(100,1);
lambda = 0.5;
act(1)=1;
for i=1:99
    act(i+1) = act(i)-lambda*(act(i)-exc);
end
% exc=0;
% for i=61:70
%     act(i+1) = act(i)-lambda*(act(i)-exc);
% end
% exc=0.3;
% for i=71:99
%     act(i+1) = act(i)-lambda*(act(i)-exc);
% end
figure; plot(act)