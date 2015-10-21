% $Id: loadOBJ.m 767 2012-04-26 20:19:41Z aaron $
% $URL: svn+ssh://svn.tau-tek.com/home/svn/FaST/trunk/Models/loadOBJ.m $
%
% Load an obj model into a structure
function a=loadOBJ(fNam)
[fd,msg] = fopen(fNam,'rt');
a=struct();
groupInfo.g="UNINITIALIZED";
groupInfo.usemtl = "UNINITIALIZED";
groupInfo.SurfaceCode = -99;
groupInfo.BulkCode=-99;
groupInfo.Temp=-99;
groupInfo.fIndex = 0;
prevGroupInfo = groupInfo;
a.GroupInfo = groupInfo;
if (fd == -1)
  disp(msg);
  return;
end
nuv=0;
nv=0;
nvn=0;
nvt=0;
nf=0;
lineNo=0
while(1)
   line = fgets(fd);
   lineNo = lineNo+1;
   if (~ischar(line) && (line < 0))
      return;
   end
%disp(sprintf('line %d\t:%s:\n',lineNo,line));
   groupInfo = parseGroupInfo(line,groupInfo);
   line = compactLine(line);
if (mod(lineNo,1000) == 0), fprintf(1,'%d)%s(\r',lineNo,line); end
   len = length(line);
   if (len > 5)
      if (strcmp(line(1:2),'v '))
         v = sscanf(line(3:len),'%f');
         nv=nv+1;
         a.v(nv,:) = v';
      elseif (strcmp(line(1:3),'vt '))
         vt = sscanf(line(4:len),'%f');
         nvt=nvt+1;
         a.vt(nvt,:) = vt';
      elseif (strcmp(line(1:3),'vn '))
         vn = sscanf(line(4:len),'%f');
         nvn=nvn+1;
         a.vn(nvn,:) = vn';
      elseif (strcmp(line(1:2),'o '))
         nuv=nuv+1;
         disp(sprintf('Texture %d at line %d : %s', nuv, lineNo, line));
      elseif (strcmp(line(1:2),'f '))
         if (~groupInfoMatch(groupInfo,prevGroupInfo))
            groupInfo.fIndex = nf;
            prevGroupInfo = groupInfo;
            a.GroupInfo(length(a.GroupInfo)+1)=groupInfo;
            groupInfo
            disp(sprintf('at line %d',lineNo));
         end
         nf=nf+1;
         %f = sscanf(line(3:len),'%d/%d/%d %d/%d/%d %d/%d/%d');
         %a.f(nf).v  = f([1 4 7])';
         %a.f(nf).vt = f([2 5 8])';
         %a.f(nf).vn = f([3 6 9])';
         [ifv,ifvt,ifvn] = parseFaces(line(3:len));
%keyboard
         a.f(nf).v  = ifv;
         a.f(nf).vt = ifvt;
         a.f(nf).vn = ifvn;
         a.f(nf).uv=int16(nuv);
      else
         disp(sprintf('%d Skipped : %s',lineNo,line));
      end
   end
end
end

function [iv,ivt,ivn] = parseFaces(lineText)
  [fa,fb,fc,count,errMsg] = sscanf(lineText,'%s %s %s','C');
  ifv=[0 0 0]; ifvt=ifv; ifvn=ifv;
  if (count ~= 3)
    disp("ERROR : All faces must be triangles.");
    return;
  end

  [v,vt,vn] = parseFaceVertexIndices(fa);  iv(1)=v; ivt(1)=vt; ivn(1)=vn;
  [v,vt,vn] = parseFaceVertexIndices(fb);  iv(2)=v; ivt(2)=vt; ivn(2)=vn;
  [v,vt,vn] = parseFaceVertexIndices(fc);  iv(3)=v; ivt(3)=vt; ivn(3)=vn;
end

function [v,vt,vn] = parseFaceVertexIndices(fs)
  v=0; vt=0; vn=0;
  sp = findstr(fs,'/');
  if (length(sp) < 1)
    v = str2num(fs);
    return;
  end
  if (length(sp) < 2)
    ss = fs(1:sp(1)-1);    v  = str2num(ss);
    if (length(fs) >= sp(1)+1)
      ss = fs(sp(1)+1:end);  vt = str2num(ss);
    end
    return;
  end

  ss = fs(1:sp(1)-1);    v  = str2num(ss);
  if (sp(2) > sp(1)+1)
    ss = fs(sp(1)+1:sp(2)-1);  vt = str2num(ss);
  end
  if (sp(2)+1 <= length(fs))
    ss = fs(sp(2)+1:end);  vn = str2num(ss);
  end

end

function b=compactLine(aa)
a=aa;
b='';
i = index(a,'#');
if (i > 0)
   if (i == 1)
      return;
   end
   a = a(1:i-1);  % strip comments
end
aLen=length(a);
while((aLen > 0) && isspace(a(aLen)) )
   a = a(1:aLen-1);  % trim trailing spaces
   aLen=aLen-1;
end
n=length(a);
ib=0;
i=1;
while(i <= n)
   while((i <= n) && isspace(a(i)) )
      i=i+1; % trim leading space
      if (i>n), return; end
   end
   while((i <= n) && ~isspace(a(i)) )
      ib=ib+1;
      b(ib) = a(i);
      i=i+1;
      if (i>n), return; end
   end
   if (i < n)
      ib=ib+1;
      b(ib)=' '; % add word delimiter
   end
   i=i+1;
end
end

% check if group information is in this line, if so, update current group info struct
function gi = parseGroupInfo(line,gi0)
gi=gi0;
key = 'g '; len=length(key);
if (strncmp(key,line,len))
   gi.g = compactLine(line(len+1:length(line)));
end
key = 'usemtl '; len=length(key);
if (strncmp(key,line,len))
   gi.usemtl = compactLine(line(len+1:length(line)));
end
key = '#Surface Code: '; len=length(key);
if (strncmp(key,line,len))
   gi.SurfaceCode = str2num(compactLine(line(len+1:length(line))));
end
key = '#Bulk Code: '; len=length(key);
if (strncmp(key,line,len))
   gi.BulkCode = str2num(compactLine(line(len+1:length(line))));
end
key = '#Initial Temperature (Degrees Kelvin): ';  len=length(key);
if (strncmp(key,line,len))
   gi.Temp = str2num(compactLine(line(len+1:length(line))));
end
end

% return true if these two groups are essentially the same
function b = groupInfoMatch(g1,g2)
b=1; % assume match, they usually will, unless there was a significant change
if (~strcmp(g1.g,g2.g)), b=0; end
if (~strcmp(g1.usemtl,g2.usemtl)), b=0; end
if (abs(g1.SurfaceCode-g2.SurfaceCode)>.1), b=0; end
if (abs(g1.BulkCode-g2.BulkCode)>.1), b=0; end
if (abs(g1.Temp-g2.Temp)>1), b=0; end
end
