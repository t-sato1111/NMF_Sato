clear; close all; clc;

[audioSig, fs] = audioread("pfcl.wav"); % 音声データ取り出し
sigLen = length(audioSig); % 音声データの長さ
F = DGTtool('windowShift',256,'windowLength',1024,'windowName','hann'); % DGTtoolオブジェクトを作成
X = abs(F(audioSig)); % STFTしたデータの位相情報を取り除いて保存
theta = angle(F(audioSig)); % STFTしたデータの位相情報を保存

K = 2; % 基底の数
doCount = 200; % 反復回数

[I, J] = size(X); % Xの行の数をI、列の数をJに代入
W = rand(I, K); % 各要素がランダムな正の整数のI×K行列
H = rand(K, J); % 各要素がランダムな正の整数のK×J行列
Loss = zeros(doCount, 1); % 類似度の計算に使用

for l = 0 : doCount - 1 % WとHを反復更新
    W = W .* (X * H.') ./ (W * (H * H.')); % Wの更新式
    W = max(W, eps); % 0割りを阻止
    H = H .* (W.' * X) ./ ((W.' * W) * H); % Hの更新式
    H = max(H, eps); % 0割りを阻止

    Y = W * H; % Xを復元
    Loss(l + 1, 1) = sum((abs(X - Y)) .^ 2, "all"); % 二乗Euclid距離を計算
end

semilogy(Loss); % 縦軸は対数で二乗Euclid距離を表示
mixSig = F.pinv(Y .* exp(1i * theta)); % 復元した行列の各要素に位相情報を追加
sig = zeros(sigLen, K); % 分離した要素を入れるための行列

for i = 1:K % 基底の数だけ繰り返し
    tmp = F.pinv((W(:, i) * H(i, :)) .* exp(1i * theta)); % 要素を計算して位相情報を追加
    sig(:, i) = tmp(1:sigLen, :); % 計算した要素を、ゼロパディングで長くなった部分を取り除いてから保存
end

% 音の出力
outputDir = "./output/"; % 出力先を指定
if ~exist(outputDir, 'dir') %指定のフォルダがなければ作る
    mkdir(outputDir);
end

audiowrite(outputDir+"mixSig.wav", mixSig, fs); % iSTFTした音声ファイル
for i = 1:K % 基底の数だけ繰り返し
audiowrite(outputDir+"sig"+ i + ".wav", sig(:, i), fs); % 分離した各要素の音声ファイル
end