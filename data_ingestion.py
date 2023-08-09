import argparse
import logging
function q = computweight(particle, z_RD, R_Resolution, V_Resolution, Xx, Yy, delta_n,Nr,Nd)
q = 1;
p=1;
% particle是粒子的状态，包含了距离和速度两个部分
range_particle = particle(1);
velocity_particle = particle(2);

% 根据粒子的状态，计算对应的距离和速度区间的索引
range_index = ceil((range_particle - min(Xx)) / R_Resolution);
velocity_index = ceil((velocity_particle - min(Yy)) / V_Resolution);

% 检查索引是否有效，如果不在有效的范围内，则给粒子一个最小的权重
if range_index < 1 || range_index > length(Xx) || velocity_index < 1 || velocity_index > length(Yy)
    q = 1; % 无效索引对应的最小权重
    return;
end

% 从z_RD中提取出与粒子状态对应的雷达测量值
measurement = z_RD(range_index, velocity_index);
% % 计算R的似然性
% R_likelihood = exp(-((range_index+214)*R_Resolution - range_particle)^2 / (2*delta_n^2)) / sqrt(2*pi*delta_n^2);
% % 计算V的似然性
% V_likelihood = exp(-((velocity_index-19)*V_Resolution - velocity_particle)^2 / (2*delta_n^2)) / sqrt(2*pi*delta_n^2);
% h = R_likelihood * V_likelihood;
%
% % 粒子的权重与幅度似然性和相位似然性的乘积成正比
h = 1;
q = q*(2*delta_n^2/h)*exp((0.5/delta_n^2-1/h)*abs(measurement)^2);
h9=[];
for ii=max(1,range_index-p):min(Nr,range_index+p)
    for jj=max(1,velocity_index-p):min(Nd,velocity_index+p)
%         根据高斯噪声模型计算距离的似然性
        R_likelihood = exp(-((ii+213)*R_Resolution - range_particle)^2 / (2*delta_n^2));
%         计算速度的似然性
        V_likelihood = exp(-((jj-22)*V_Resolution - velocity_particle)^2 / (2*delta_n^2));
        h = R_likelihood * V_likelihood;
        h9=[h9 h];
%         h_rdb_p=Amp^2*exp(-Lr*(r(ii)-zp(1))^2/Dr-Ld*(d(jj)-zp(2))^2/Dd-Lb*(b(mm)-zp(3))^2/Db)+2*delta_n^2;
        % h = 1;
        % q=q*(2*delta_n^2/h)*exp((0.5/delta_n^2-1/h)*abs(measurement)^2);
        q=q*(2*delta_n^2/(h+2*delta_n^2))*exp((0.5/delta_n^2-1/(h+2*delta_n^2))*abs(measurement)^2);
    end
end

end
function q = computweight(particle, z_RD, R_Resolution, V_Resolution, Xx, Yy, delta_v,Nr,Nd)
q = 1
p = 1
% particle是粒子的状态，包含了距离和速度两个部分
range_particle = particle(1);
velocity_particle = particle(2);

% 根据粒子的状态，计算对应的距离和速度区间的索引
range_index = ceil((range_particle - min(Xx)) / R_Resolution);
velocity_index = ceil((velocity_particle - min(Yy)) / V_Resolution);
% 从z_RD中提取出与粒子状态对应的雷达测量值

% 检查索引是否有效，如果不在有效的范围内，则给粒子一个最小的权重
if range_index < 1 || range_index > length(Xx) || velocity_index < 1 || velocity_index > length(Yy)
    q = 1e-6; % 无效索引对应的最小权重
    return;
end
measurement = z_RD(range_index, velocity_index);
for ii=max(1,range_index-p):min(Nr,range_index+p)
    for jj=max(1,velocity_index-p):min(Nd,velocity_index+p)
        % 根据高斯噪声模型计算距离的似然性
        R_likelihood = exp(-((ii+214)*R_Resolution - range_particle)^2 / (2*delta_v^2)) / sqrt(2*pi*delta_v^2);
        % 计算速度的似然性
        V_likelihood = exp(-((jj-19)*V_Resolution - velocity_particle)^2 / (2*delta_v^2)) / sqrt(2*pi*delta_v^2);
        h = R_likelihood * V_likelihood;
        % h_rdb_p=Amp^2*exp(-Lr*(r(ii)-zp(1))^2/Dr-Ld*(d(jj)-zp(2))^2/Dd-Lb*(b(mm)-zp(3))^2/Db)+2*delta_n^2;
        q=q*(2*delta_v^2/h)*exp((0.5/delta_v^2-1/h)*abs(measurement)^2);
    end
end

end
from autogpt.commands.file_operations import ingest_file, search_files
from autogpt.config import Config
from autogpt.memory import get_memory

cfg = Config()


def configure_logging():
    logging.basicConfig(
        filename="log-ingestion.txt",
        filemode="a",
        format="%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        level=logging.DEBUG,
    )
    return logging.getLogger("AutoGPT-Ingestion")


def ingest_directory(directory, memory, args):
    """
    Ingest all files in a directory by calling the ingest_file function for each file.

    :param directory: The directory containing the files to ingest
    :param memory: An object with an add() method to store the chunks in memory
    """
    try:
        files = search_files(directory)
        for file in files:
            ingest_file(file, memory, args.max_length, args.overlap)
    except Exception as e:
        print(f"Error while ingesting directory '{directory}': {str(e)}")


def main() -> None:
    logger = configure_logging()

    parser = argparse.ArgumentParser(
        description="Ingest a file or a directory with multiple files into memory. "
        "Make sure to set your .env before running this script."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--file", type=str, help="The file to ingest.")
    group.add_argument(
        "--dir", type=str, help="The directory containing the files to ingest."
    )
    parser.add_argument(
        "--init",
        action="store_true",
        help="Init the memory and wipe its content (default: False)",
        default=False,
    )
    parser.add_argument(
        "--overlap",
        type=int,
        help="The overlap size between chunks when ingesting files (default: 200)",
        default=200,
    )
    parser.add_argument(
        "--max_length",
        type=int,
        help="The max_length of each chunk when ingesting files (default: 4000)",
        default=4000,
    )

    args = parser.parse_args()

    # Initialize memory
    memory = get_memory(cfg, init=args.init)
    print("Using memory of type: " + memory.__class__.__name__)

    if args.file:
        try:
            ingest_file(args.file, memory, args.max_length, args.overlap)
            print(f"File '{args.file}' ingested successfully.")
        except Exception as e:
            logger.error(f"Error while ingesting file '{args.file}': {str(e)}")
            print(f"Error while ingesting file '{args.file}': {str(e)}")
    elif args.dir:
        try:
            ingest_directory(args.dir, memory, args)
            print(f"Directory '{args.dir}' ingested successfully.")
        except Exception as e:
            logger.error(f"Error while ingesting directory '{args.dir}': {str(e)}")
            print(f"Error while ingesting directory '{args.dir}': {str(e)}")
    else:
        print(
            "Please provide either a file path (--file) or a directory name (--dir)"
            " inside the auto_gpt_workspace directory as input."
        )


if __name__ == "__main__":
    main()
clc;
clear;
load RDM.mat;
R_Trx = 0;
R_rangeMin = 520+R_Trx;
R_rangeMax = 580+R_Trx;
R_Resolution = 2.44;
R_binMin = floor(R_rangeMin/R_Resolution)+1;
R_binMax = floor(R_rangeMax/R_Resolution)+1;
Nr = R_binMax - R_binMin;
RbinRange = R_binMin:R_binMax;
Xx = (RbinRange-1)*R_Resolution-R_Trx;   %meter

PRI          = 2.5e-3;   %雷达脉冲重复间隔
fc           = 3.7e9;
frame_length     =128;   %转换后的帧长
v_range  = -19:20; 
V_Resolution = 3e8/fc/PRI/frame_length*3.6;
V_rangeMin = v_range(1); 
V_rangeMax = v_range(end);  
V_binMin = floor(V_rangeMin/V_Resolution) + frame_length/2;
V_binMax = floor(V_rangeMax/V_Resolution) + frame_length/2;
Nd = V_binMax - V_binMin;
VbinRange = V_binMin:V_binMax;
Yy = (VbinRange-1-frame_length/2)*V_Resolution;    %kmph
%% parameters initialization
steps = 30;
x = zeros(2,steps); %状态
z = zeros(2,steps); %测量
T = 0.1;
F_cv=[1 -T ; % transition matrix of cv model
      0 1];
G=[T^2/2 0; % process noise transition matrix, assuming the noise distribution is the same for x and y axis
    T 0];
delta_v=0.01; %process noise standard deviation here, equal to 0;
Qv=G*delta_v^2*G';% process noise covariance
h_RD = zeros(Nr,Nd,steps); % signal contribution to the cell
z_rdb_I=zeros(Nr,Nd,steps); % real part of z (complex amplitude)
z_rdb_Q=zeros(Nr,Nd,steps); % imaginary part of z
z_RD = RDM(RbinRange,VbinRange,1:steps); % ROI量测值
fai=0;
N=1000; % number of particles
xp=zeros(2,steps,N); % state particle
zp=zeros(2,steps,N); % range,doppler,bearing particl
xe=zeros(2,steps); % state estimate
axe=zeros(2,steps);
Pb=zeros(1,steps); % existing probability
E=zeros(steps,N); % existing variable
q=ones(steps,N); % weight 
EE=[0.9 0.1 % E11,E10,E01,E00
    0.1 0.9];
delta_v_p=1.5;

% particles initialization
for i = 1:N
    xp(1,1,i)=R_rangeMin+(R_rangeMax-R_rangeMin)*rand;
    xp(2,1,i)=V_rangeMin+(V_rangeMax-V_rangeMin)*rand;
    q(1,i) = 1/N;
    if rand<0.1
        E(1,i)=1;
    else
        E(1,i)=0;
    end
    xe(:,1)=xe(:,1)+E(1,i)*xp(:,1,i); % initial estimation
%     axe(:,1)=axe(:,1)+E(1,i)*xp(:,1,i);
    Pb(1,1)=Pb(1,1)+E(1,i); % initial probability
end
Pb(1,1)=Pb(1,1)/N;
axp=xp;
%% iteration starts
steps_temp=steps; %迭代开始
for k = 1 : steps_temp - 1;
    for i = 1:N
        if E(k,i)==0    
            if rand<0.1
                E(k+1,i)=1;
            else
                E(k+1,i)=0;
            end
        else
            if rand<0.1
                E(k+1,i)=0;
            else
                E(k+1,i)=1;
            end
        end
         if E(k,i)==1&&E(k+1,i)==1   % paricle maintain 
            axp(:,k+1,i)=F_cv*axp(:,k,i)+G*delta_v_p*randn(2,1);    
        elseif E(k,i)==0&&E(k+1,i)==1 % new particle, note: we can choose a better zone for new particles 
            axp(1,k+1,i)=R_rangeMin+(R_rangeMax-R_rangeMin)*rand;
            axp(2,k+1,i)=V_rangeMin+(V_rangeMax-V_rangeMin)*rand;
        elseif E(k+1,i)==0
            axp(:,k+1,i)=0;
         end
         xp(:,k+1,i)=axp(:,k+1,i);
         q(k+1,i) = 1/N; % TODO compute weights
    end
    qsum = sum(q(k+1,:));
    q(k+1,:) = q(k+1,:)/qsum;
    
    % resample
    xp_temp=xp;
    axp_temp=axp;
    EE_temp=E(k+1,:);
    for i=1:N
        uu=rand;
        qq_sum=0;
        for j=1:N
            qq_sum=qq_sum+q(k+1,j);
            if qq_sum>=uu
               xp(:,k+1,i)=xp_temp(:,k+1,j);
                axp(:,k+1,i)=axp_temp(:,k+1,j);
                E(k+1,i)=EE_temp(j);
                break;
            end
        end
    end
    for i=1:N
        Pb(1,k+1)=Pb(1,k+1)+E(k+1,i);
%          xe(:,k+1)=xe(:,k+1)+E(k+1,i)*xp(:,k+1,i);
         axe(:,k+1)=axe(:,k+1)+E(k+1,i)*axp(:,k+1,i);
    end
    if Pb(1,k+1)~=0;
%        xe(:,k+1)=xe(:,k+1)/Pb(1,k+1);
       axe(:,k+1)=axe(:,k+1)/Pb(1,k+1);%估计状态
    end

    Pb(1,k+1)=Pb(1,k+1)/N;
end % pf迭代结束

% 记录目标出现概率与位置
PT=zeros(1,steps);
tt=0;
for k=1:steps
    if Pb(k)>0.6;
        PT(k)=1;
        tt=tt+1;
        Fxe(:,tt)=axe(:,k);
    end
end

% 画图
figure(1);
plot(1:steps_temp,Pb(1:steps_temp),'ro-');
title('Probability');
