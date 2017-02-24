clear
global coeff;

for f=1:200


    r=.33;
    gamma=.011;
    koff=.0045;
    kon=.00027*10^(f/40);

    XD(1) = 1; %XD  --> 0 or 1; I or A 
    XR(1) = 1; %XR  --> mRNA # 

    c=[koff kon r gamma];     %R1: I -> A, R2: A -> I, R3: m -> m+1, R4: m -> m-1 

    MaxT=1000;
    for i=1:MaxT

        h=[1-XD(i) XD(i) XD(i) XR(i)];
        a=h.*c;
        a0=sum(a); a1=a(1); a2=a(1)+a(2); a3=a(1)+a(2)+a(3);

        Random1 = rand;
        Random2 = rand;

        tau(i) = (1/a0)*log(1/Random1);

        if 0 < Random2*a0 && Random2*a0 <= a1
            XD(i+1) = 1-XD(i);
            XR(i+1) = XR(i);
        elseif a1 < Random2*a0 && Random2*a0 <= (a2)
            XD(i+1) = 1-XD(i);
            XR(i+1) = XR(i);
        elseif (a2) < Random2*a0 && Random2*a0 <= (a3)

            XD(i+1) = XD(i);
            XR(i+1) = XR(i)+1;
        elseif (a3) < Random2*a0 && Random2*a0<= a0
            XD(i+1) = XD(i);
            XR(i+1) = XR(i)-1;
        end
    end



    MeanM=sum(XR(1:end-1).*tau)/sum(tau);
    SecondMoment=sum((XR(1:end-1).^2).*tau)/sum(tau);

    VarianceM=SecondMoment-MeanM^2;
    Fanofactor(f)=(SecondMoment-(MeanM.^2))/MeanM;

    fcin(f) = (1+kon/koff)^(-1);
    kof(f) = koff;
    konn(f)= kon;
    Mean_M(f) = MeanM;
    Mean_M_Max = max(Mean_M(f));

end 


scatter(fcin,Fanofactor);
set(gca,'xscale','log');

Curve_fit_function(fcin, Fanofactor)
list_of_coeff = coeff;


%list_of_coeff
koff/gamma
%Error_percentage = 100 * (list_of_coeff(1:end-1)-(koff/gamma))/(koff/gamma);