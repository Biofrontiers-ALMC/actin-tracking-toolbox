
%reader = BioformatsImage('WTbeta_0,5spf_1.nd2');

Isum = zeros(reader.height, reader.width);

for iT = 1:reader.sizeT
    Isum = mean(cat(3, Isum, double(getPlane(reader, 1, 1, iT))), 3);
end
%%
%Take a rough mask to identify objects
Isum = Isum ./ max(Isum(:));

mask = imbinarize(Isum, 'adaptive');

%Use the detected pixels to fit a 2D gaussian
[row, col] = find(mask);

%Fit measured data to a 2D Gaussian
gauss2D =  fittype('A * exp(-((xx - B).^2 + (yy - C).^2) / (2*D.^2))',...
    'independent', {'xx', 'yy'});

%xyFit = cat(1, rp.Centroid);
xyFit = [col, row];
intFit = cat(1, Isum(mask));

[maxInt, maxIntInd] = max(intFit);

initGuess = [maxInt, xyFit(maxIntInd, 1), xyFit(maxIntInd, 2), 1000];

fitObj = fit(xyFit, intFit, gauss2D, 'StartPoint', initGuess);

figure;
plot(fitObj, xyFit, intFit)

%Reconstruct the intensity profile
xx = 1:size(Isum, 2);
yy = 1:size(Isum, 1);
[xx, yy] = meshgrid(xx, yy);

reconInt = fitObj.A * exp(-((xx - fitObj.B).^2 + (yy - fitObj.C).^2) / (2*fitObj.D.^2));

figure;
imshow(reconInt, [])

Itest = getPlane(reader, 1, 1, 1);

Icorr = double(Itest) ./ reconInt;
imshowpair(Itest, Icorr, 'montage')
