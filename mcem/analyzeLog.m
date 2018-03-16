function res = analyseLog(mcLogFile, icmLogFile, destDir)
% analyseLog(mcLogFile, icmLogFile)

[icm.tag, icm.emIter, icm.beta, icm.priorll, icm.cll, icm.jointll, icm.mu1, icm.sigma1, icm.numPoints1, icm.mu2, icm.sigma2, icm.numPoints2, icm.mu3, icm.sigma3, icm.numPoints3, icm.mu4, icm.sigma4, icm.numPoints4] = textread(icmLogFile, '%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'headerlines', 0);
[mc.tag, mc.emIter, mc.beta, mc.priorll, mc.cll, mc.jointll, mc.mu1, mc.sigma1, mc.numPoints1, mc.mu2, mc.sigma2, mc.numPoints2, mc.mu3, mc.sigma3, mc.numPoints3, mc.mu4, mc.sigma4, mc.numPoints4] = textread(mcLogFile, '%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'headerlines', 0);

if ~exist(destDir, 'dir')
    mkdir(destDir);
end;
    
figure(1);
plot(icm.emIter, icm.beta, 'r*-', mc.emIter, mc.beta, 'b*-', 'LineWidth', 3);
xlabel('EM iteration');
ylabel('\beta');
%saveas(1, [ destDir, '/'comp_beta.png');
saveas(1, [ destDir, '/', 'comp_beta.png']);
 
figure(1);
plot(icm.emIter, icm.jointll, 'r*-', mc.emIter, mc.jointll, 'b*-', 'LineWidth', 3);
xlabel('EM iteration');
ylabel('Joint log-likelihood');
saveas(1, [ destDir, '/', 'comp_jointll.png']);

figure(1);
plot(icm.emIter, icm.mu1, 'r*-', mc.emIter, mc.mu1, 'b*-', 'LineWidth', 3);
xlabel('EM iteration');
ylabel('mu1');
saveas(1, [ destDir, '/', 'comp_mu1.png']);

figure(1);
plot(icm.emIter, icm.mu2, 'r*-', mc.emIter, mc.mu2, 'b*-', 'LineWidth', 3);
xlabel('EM iteration');
ylabel('mu2');
saveas(1, [ destDir, '/', 'comp_mu2.png']);

figure(1);
plot(icm.emIter, icm.mu3, 'r*-', mc.emIter, mc.mu3, 'b*-', 'LineWidth', 3);
xlabel('EM iteration');
ylabel('mu3');
saveas(1, [ destDir, '/', 'comp_mu3.png']);

figure(1);
plot(icm.emIter, icm.mu4, 'r*-', mc.emIter, mc.mu4, 'b*-', 'LineWidth', 3);
xlabel('EM iteration');
ylabel('mu4');
saveas(1, [ destDir, '/', 'comp_mu4.png']);

figure(1);
plot(icm.emIter, icm.sigma1, 'r*-', mc.emIter, mc.sigma1, 'b*-', 'LineWidth', 3);
xlabel('EM iteration');
ylabel('sigma1');
saveas(1, [ destDir, '/', 'comp_sigma1.png']);