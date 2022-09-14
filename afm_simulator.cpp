#include "points3d.hpp"
#include <cmath>


class afm{
  private:
  points3d surface;
  points3d probe;
  float max_surface = -1000000;
  float min_probe = 1000000;
  float DX = 0;
  float DY = 0;
  float LennardJonesForce(vector<vector<float>> p, vector<vector<float>> s, float dx, float dy, float z);
  public:
  afm(points3d s, points3d p) :surface(s),probe(p){
    // probeの下限の取得
    for(auto &p_point : probe.points){
      if(min_probe > p_point[2]){
         min_probe = p_point[2];
         DX = p_point[0];
         DY = p_point[1];
      }
    }
    // surface の上限の取得
    for(auto &s_point : surface.points){
      if(max_surface < s_point[2]){
        max_surface = s_point[2];
      }
    }
  };
  // 画像出力系
  //void scanXY(float z);
  void scanXYZ();
  void FMscanXYZ();
  void scanXY(float z);
  void scanZ(float x, float y);
};

// 制御なしのconstant hightでの計測
void afm::scanXY(float z){
  int nx = 25;
  int ny = 25;
  int nz = 30;
  float dx = 1;
  float dy = 1;
  float dz = 0.1;
  cout<<surface.points.size()<<" "<<probe.points.size()<<endl;
  cout<<nx<<" "<<ny<<" "<<nz<<endl;
  cout<<dx<<" "<<dy<<" "<<dz<<endl;
  for(auto point : surface.points){
    cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<endl;
  }
  for(auto point : probe.points){
    cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<endl;
  }
  vector<vector<float>> img(nx,vector<float>(ny));
  for(int x=0;x<nx;x++){
    for(int y=0;y<ny;y++){
      float f = LennardJonesForce(probe.points, surface.points, x*dx-DX, y*dy-DY, z);
      cout<<x<<" "<<y<<" "<<z<<" "<<f<<endl;
    }
  }
};

void afm::scanZ(float x, float y){
  float dz = 0.1;
  int z_s = 30;
  int z_e = 80;
  cout<<50<<endl;
  for(int i=z_s;i<z_e;i++){
    float f = LennardJonesForce(probe.points, surface.points, x-DX, y-DY, i*dz);
    cout<<f<<endl;
  }
};

//表面から1~6の範囲でスキャン
void afm::scanXYZ(){
  int nx = 50;
  int ny = 50;
  int nz = 30;
  float dx = 1;
  float dy = 1;
  float dz = 0.1;
  //surface原子の位置;
  cout<<surface.points.size()<<" "<<probe.points.size()<<endl;
  cout<<nx<<" "<<ny<<" "<<nz<<endl;
  cout<<dx<<" "<<dy<<" "<<dz<<endl;
  for(auto point : surface.points){
    cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<endl;
  }
  for(auto point : probe.points){
    cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<endl;
  }
  for(int z=20;z<nz+20;z++){
    for(int y=0;y<ny;y++){
      for(int x=0;x<nx;x++){
        float f = LennardJonesForce(probe.points, surface.points, x*dx-DX, y*dy-DY, z*dz);
        cout<<x<<" "<<y<<" "<<z<<" "<<f<<endl;
      }
    }
  }
};

void afm::FMscanXYZ(){
  int nx = 100;
  int ny = 100;
  int nz = 30;
  float dx = 0.5;
  float dy = 0.5;
  float dz = 0.1;
  float dA = 0.01;
  //surface原子の位置;
  cout<<surface.points.size()<<" "<<probe.points.size()<<endl;
  cout<<nx<<" "<<ny<<" "<<nz<<endl;
  cout<<dx<<" "<<dy<<" "<<dz<<endl;
  for(auto point : surface.points){
    cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<endl;
  }
  for(auto point : probe.points){
    cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<endl;
  }
  for(int z=20;z<nz+20;z++){
    for(int y=0;y<ny;y++){
      for(int x=0;x<nx;x++){
        float f1 = LennardJonesForce(probe.points, surface.points, x*dx-DX, y*dy-DY, z*dz+dA/2.0);
        float f2 = LennardJonesForce(probe.points, surface.points, x*dx-DX, y*dy-DY, z*dz-dA/2.0);
        cout<<x<<" "<<y<<" "<<z<<" "<<f1/dA-f2/dA<<endl;
      }
    }
  }
};

float afm::LennardJonesForce(vector<vector<float>> probe, vector<vector<float>> surface, float dx, float dy, float z){
  float F = 0;
  float eV = 1.609e-19;
  float k = 1.38062e-23;
  float m0 = 1e-26;
  float d = 1e-10;
  float sigma = 2.951;
  float epsilon = 0.2294;
  for(auto p: probe){
    // プローブ先の移動を擬似的に表現
    p[0] += dx;
    p[1] += dy;
    // 高さ制御
    p[2] = p[2] - min_probe + z + max_surface;
    for(auto s: surface){
      float distance = pow(pow(p[0] - s[0], 2) + pow(p[1] - s[1], 2) + pow(p[2] - s[2], 2), 0.5);// 純粋な距離
      float distance_z = p[2] - s[2];
      // LennardJonesの式に入れて計算
      //F = -48*epsilon*((sigma/r)**12 - 0.5 * (sigma/r)**6)*R[2]/(r**2)
      F += -48 * epsilon * (pow(sigma / distance,12) - 0.5 * pow(sigma / distance, 6)) * distance_z / pow(distance, 2);
    }
  }
  return F;
}


// 曲率
int main(int argc, char *argv[]){
  points3d surface(0, false, 4.078);
  points3d probe(0, false, 4.078);

  surface.add_atom(11,14,0);
  surface.add_atom(11,11,0);
  surface.add_atom(14,14,0);
  surface.add_atom(14,11,0);
  string option = argv[1];
  if(option=="xy"){
    const float probe_r = stof(argv[3]); //曲率
    if(probe_r==0){
      probe.add_atom(0,0,0);
    }else{
      probe.make_probe(probe_r);
      int needle_num;cin>>needle_num;
      for(int i=0;i<needle_num;i++){
        float x,y,z;
        cin>>x>>y>>z;
        probe.add_atom(x,y,z);
      }
    }
    afm AFM(surface, probe);
    AFM.scanXY(stof(argv[2]));
  }else if(option=="z"){
    const float probe_r = stof(argv[4]); //曲率
    if(probe_r==0){
      probe.add_atom(0,0,0);
    }else{
      probe.make_probe(probe_r);
      int needle_num;cin>>needle_num;
      for(int i=0;i<needle_num;i++){
        float x,y,z;
        cin>>x>>y>>z;
        probe.add_atom(x,y,z);
      }
    }
    afm AFM(surface, probe);
    AFM.scanZ(stof(argv[2]),stof(argv[3]));
  }
}
//三次元スキャン
//>more needle\100_15_5_1.txt | a.exe 100|py run.py 100_15_5_1
//more needle\100_15_5_1.txt | a.exe 100|py run2d.py 100_15_5_1