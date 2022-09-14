#ifndef POINTS3D
#define POINTS3D
#include <vector>
#include <stdio.h>
#include <random>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <functional>
#include <algorithm>
using namespace std;
class points3d{
  private:

  public:
  // 原子の種類,乱数を使うか,ユニットセルの長さ
  int atom_kind;
  bool random_atom;
  float unit_cell_length;
  vector<vector<float>> points;

  // 原子の種類、乱数を使うかの設定
  points3d(int ak, bool ra, float ucl): atom_kind(ak),random_atom(ra),unit_cell_length(ucl) {};
  // 平面を生成する
  void face_centered_cubic(float x, float y, float z);
  void print();
  void shift(float x, float y, float z);
  void trim(function<bool(tuple<float,float,float>)> fn);
  void make_probe(float r);
  void add_atom(float x, float y, float z);
};

//原子を追加する
void points3d::add_atom(float x, float y, float z){
  points.push_back({x, y, z});
}

// [0A,xA]の範囲で埋めて行く
void points3d::face_centered_cubic(float x, float y, float z){
  int lx = x / unit_cell_length;
  int ly = y / unit_cell_length;
  int lz = z / unit_cell_length;
  for(int ix = 0; ix<= lx; ix++){
    for(int iy = 0; iy<=ly; iy++){
      for(int iz = 0; iz<=lz; iz++){
        if(random_atom){
          random_device rd;
          default_random_engine eng(rd());
          uniform_real_distribution<float> distr(-0.1, 0.1);
          points.push_back({unit_cell_length * ix + distr(eng), unit_cell_length * iy + distr(eng), unit_cell_length * iz+ distr(eng)});
          points.push_back({unit_cell_length * (float)(ix + 0.5)+ distr(eng), unit_cell_length * (float)(iy + 0.5) + distr(eng), unit_cell_length * iz+ distr(eng)});
          points.push_back({unit_cell_length * (float)(ix + 0.5)+ distr(eng), unit_cell_length * iy+ distr(eng) , unit_cell_length * (float)(iz + 0.5)+ distr(eng)});
          points.push_back({unit_cell_length * ix+ distr(eng), unit_cell_length * (float)(iy + 0.5)+ distr(eng) , unit_cell_length * (float)(iz + 0.5)+ distr(eng)});
        }else{
          points.push_back({unit_cell_length * ix, unit_cell_length * iy , unit_cell_length * iz});
          points.push_back({unit_cell_length * (float)(ix + 0.5), unit_cell_length * (float)(iy + 0.5) , unit_cell_length * iz});
          points.push_back({unit_cell_length * (float)(ix + 0.5), unit_cell_length * iy , unit_cell_length * (float)(iz + 0.5)});
          points.push_back({unit_cell_length * ix, unit_cell_length * (float)(iy + 0.5) , unit_cell_length * (float)(iz + 0.5)});
        }
      }
    }
  }
};

void points3d::print(){
  for(const auto&point :points){
    cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<endl;
  }
};

void points3d::shift(float x, float y, float z){
  for(auto &point : points){
    point[0] += x;
    point[1] += y;
    point[2] += z;
  }
};

void points3d::trim(function<bool(tuple<float,float,float>)> fn){
  auto itr = points.begin();
  while (itr != points.end()) {
    if(fn(tie(itr->at(0),itr->at(1),itr->at(2)))){
      ++itr;
    }else{
      itr = points.erase(itr);
    }
  }
};

void points3d::make_probe(float r){
  float under_limit = min<float>({0,15-r});
  points3d::face_centered_cubic(2*r+10,2*r+10,2*r+10); 
  points3d::shift(-r,-r,-r);
  auto itr = points.begin();
  int s = 0;
  vector<vector<float>> new_points(0);
  for(auto const point : points){
    float x = point[0];
    float y = point[1];
    float z = point[2];
    if(x*x + y*y + z*z < r*r && z<=under_limit){
      //条件に合致したものを新しいものに入れてく
      vector<float> tmp = {x,y,z};
      new_points.push_back(tmp);
    }
  }
  points = new_points;
};
#endif