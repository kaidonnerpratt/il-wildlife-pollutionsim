#include <SFML/Graphics.hpp>
#include "PerlinNoise.hpp"
#include <iostream>
#include <list>
#include <chrono>
#include <cmath>
#include <thread>
#include <vector>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <typeinfo>

using namespace std::placeholders;

const int CREATURE_MAX = 2000;
const int CREATURE_START = 1000;

const int lwidth = 500;
const int lheight = 500;

const int gdst = 50;


double grassArray[lwidth*lheight];


std::vector<int> preytokill;

template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &orig){   
    std::vector<T> ret;
    for(const auto &v: orig)
        ret.insert(ret.end(), v.begin(), v.end());                                                                                         
    return ret;
}  



double clockToMilliseconds(clock_t ticks){
    // units/(units/time) => time (seconds) * 1000 = milliseconds
    return (ticks/(double)CLOCKS_PER_SEC)*1000.0;
}

template<typename T>
class point{
    public:
        T x;
        T y;
        int id = 0;
        point(T xx, T yy){
            x=xx;
            y=yy;
        }
        point(T xx, T yy, int i){
            x=xx;
            y=yy;
            id=i;
        }
    T distance(point<T> s){
        return (T) hypot((x-s.x),(y-s.y));
    }
    double direction(point<T> s){
        return atan2((y-s.y),(x-s.x));
    }
    friend std::ostream& operator<<(std::ostream& os, const point& m) {
        os << "("<<std::to_string(m.x)<<","<<std::to_string(m.y)<<")";
        return os; // enable chaining
    }
    
};


template <int tilesize, int gridsize>
class refGrid{
    public:
        static const int tilecount = gridsize/tilesize;
        std::vector<point<double>> curgrid[tilecount][tilecount];
        refGrid(){}
    void reset(){
        for (int y = 0; y < tilecount; y++){
            for (int x = 0; x < tilecount; x++){
                curgrid[y][x] = std::vector<point<double>>();
            }
        }
    }
    template<class gst>
    void sort(gst &s){
        
        curgrid[(int)floor(s.y/tilesize)][(int)floor(s.x/tilesize)].push_back(s.toPoint());
    }
    point<double> getNearestSingle(point<double> s, int tx, int ty){
        if (tx < 0){
            return point<double>(-1,-1,-1);
        }
        if (ty < 0){
            return point<double>(-1,-1,-1);
        }
        if (tx >= tilecount){
            return point<double>(-1,-1,-1);
        }
        if (ty >= tilecount){
            return point<double>(-1,-1,-1);
        }

        std::vector<point<double>> current = curgrid[ty][tx];

        point<double> closest = point<double>(-1,-1,-1);
        double closestDist = -1;

        for (int i = 0; i < current.size(); i++){
            if (current[i].distance(s) < closestDist || closestDist == -1){
                closest = current[i];
                closestDist = current[i].distance(s);
            }
        }
        return closest;



    }
    point<double> getNearestAround(point<double> s){
        int tx = (int)floor(s.x/tilesize);
        int ty = (int)floor(s.y/tilesize);

        point<double> centclosest = getNearestSingle(s, tx, ty);
        if (centclosest.id != -1){
            return centclosest;
        }else{
            
            point<double> ac = getNearestSingle(s, tx+1, ty);
            point<double> bc = getNearestSingle(s, tx+1, ty+1);
            point<double> cc = getNearestSingle(s, tx, ty+1);
            point<double> dc = getNearestSingle(s, tx-1, ty+1);
            point<double> ec = getNearestSingle(s, tx-1, ty);
            point<double> fc = getNearestSingle(s, tx-1, ty-1);
            point<double> gc = getNearestSingle(s, tx, ty-1);
            point<double> hc = getNearestSingle(s, tx+1, ty-1);



            point<double> closest = point<double>(-1,-1,-1);
            double closestDist = -1;
            if (ac.id != -1 && (ac.distance(s) < closestDist || closestDist == -1)){
                closest = ac;
                closestDist = ac.distance(s);
            }
            if (bc.id != -1 && (bc.distance(s) < closestDist || closestDist == -1)){
                closest = bc;
                closestDist = bc.distance(s);
            }
            if (cc.id != -1 && (cc.distance(s) < closestDist || closestDist == -1)){
                closest = cc;
                closestDist = cc.distance(s);
            }
            if (dc.id != -1 && (dc.distance(s) < closestDist || closestDist == -1)){
                closest = dc;
                closestDist = dc.distance(s);
            }
            if (ec.id != -1 && (ec.distance(s) < closestDist || closestDist == -1)){
                closest = ec;
                closestDist = ec.distance(s);
            }
            if (fc.id != -1 && (fc.distance(s) < closestDist || closestDist == -1)){
                closest = fc;
                closestDist = fc.distance(s);
            }
            if (gc.id != -1 && (gc.distance(s) < closestDist || closestDist == -1)){
                closest = gc;
                closestDist = gc.distance(s);
            }
            if (hc.id != -1 && (hc.distance(s) < closestDist || closestDist == -1)){
                closest = hc;
                closestDist = hc.distance(s);
            }
            return closest;

        }

    }
    
};






template<class ct, int CREATURE_MAX>
class CreatureOrganizer : public sf::Drawable {
    public:
        
        refGrid<gdst, lwidth> gridsort = refGrid<gdst, lwidth>();
        sf::Image lakes_img;
        sf::Color lake_col;
        std::vector<int> spotsNeeded;
        std::vector<int> freeSpots;
        std::vector<std::thread> currentthreads;
        
        int pop = 0;


        ct CreatureList[CREATURE_MAX];
        CreatureOrganizer(sf::Image llakes_img, sf::Color llake_col, const int CREATURE_START){
            lakes_img = llakes_img;
            lake_col = llake_col;
            for (int i = 0; i < CREATURE_MAX; i++){
                int curx = lwidth * random()/RAND_MAX;
                int cury = lheight* random()/RAND_MAX;
                while (lakes_img.getPixel(curx, cury) == lake_col){
                    curx = lwidth * random()/RAND_MAX;
                    cury = lheight* random()/RAND_MAX;
                }
                CreatureList[i] = ct(curx, cury, i >= CREATURE_START, (double)random()/(double)RAND_MAX, i);
            }
        }
    ct& operator[] (size_t i){
        return CreatureList[i];
    }
    void tickRange(int start, int end){
        for (int i = start; i < CREATURE_MAX && i < end; i++){
            ct current = CreatureList[i];

            if(!current.isnull){
                pop+=1;
                current.tick(lakes_img,lake_col,grassArray,std::ref(spotsNeeded));
                CreatureList[i] = current;
            }
        }
    }
    void tickStep(int offset, int step, refGrid<gdst,lwidth> oppsort){
        for (int i = offset; i < CREATURE_MAX; i+=step){
            ct current = CreatureList[i];

            if(!current.isnull){
                pop+=1;
                current.tick(lakes_img,lake_col,grassArray,std::ref(spotsNeeded), oppsort);
                CreatureList[i] = current;
            }
        }
    }
    void tickThreaded(int threadcount, refGrid<gdst,lwidth> oppsort){ 
        
        freeSpots.clear();
        gridsort.reset();
        
        

        pop = 0;
        int step = ceil((double)CREATURE_MAX/(double)threadcount);
       
        for (int i = 0; i < threadcount; i++){
            currentthreads.push_back(std::thread(std::bind(&CreatureOrganizer::tickStep, this, i, threadcount, oppsort)));
        }
        
        
        for (int i = 0; i < threadcount; i++){
            if(currentthreads[i].joinable()){
                currentthreads[i].join();
            }
        }
        
        for (int i = 0; i < CREATURE_MAX; i++){
            if (!CreatureList[i].isnull){
                gridsort.sort<ct>(CreatureList[i]);
                if (CreatureList[i].cankill && CreatureList[i].wannakill != -1){
                    preytokill.push_back(CreatureList[i].wannakill);
                }
            }else{
                freeSpots.push_back(i);
            }
        }
        for (int i = 0;  !spotsNeeded.empty() && i < freeSpots.size(); i++){
            ct tmpc = CreatureList[spotsNeeded.back()];
            CreatureList[freeSpots[i]] = ct(tmpc.x, tmpc.y, false, 0, freeSpots[i]);
            spotsNeeded.pop_back();
        }
        currentthreads.clear();
        spotsNeeded.clear();
    }
    void tickAll(){
        pop=0;
        std::vector<int> old_spotsNeeded = spotsNeeded;
        spotsNeeded.clear();

        for (int i = 0; i < CREATURE_MAX; i++){
            ct current = CreatureList[i];
            if(!current.isnull){
                pop+=1;
                current.tick(lakes_img,lake_col,grassArray,std::ref(spotsNeeded));
                CreatureList[i] = current;
            }else{
                if (!old_spotsNeeded.empty()){
                    ct tmpc = CreatureList[old_spotsNeeded.back()];
                    CreatureList[i] =  ct(tmpc.x, tmpc.y, false, 0, i);
                    old_spotsNeeded.pop_back();
                }
            }
        }
    }
    private:
        void draw(sf::RenderTarget& target, sf::RenderStates state) const{
            for (int i = 0; i < CREATURE_MAX; i++){
                if (!CreatureList[i].isnull){
                    target.draw(CreatureList[i]);
                }
            }
            
        }
};

class Creature : public sf::Drawable {       
    public:
        bool cankill=false;
        bool wannakill=-1;
        int id = 0;

        double matetimer = 1;
        double age = 0;


        bool iseating = false;
        double hunger = 1;


        double x;                                               //x positon
        double y;                                               //y position
        double r = 2*M_PI*(double)random()/(double)RAND_MAX;    //rotation 

        int size = 1;

        bool isnull;        //does it exist
        Creature(){
            
        }
        Creature(double a, double b, bool c, double g){
            age=g;
            x = a;
            y = b;
            isnull = c;
        
        }
        Creature(double a, double b, bool c, double g, int i){
            age=g;
            x = a;
            y = b;
            isnull = c;
            id=i;
        }
    void tick(sf::Image lakes_img, sf::Color lake_col, double *grassArray, std::vector<int> &spotsNeeded, refGrid<gdst,lwidth> oppsort){

    }
    void eat(refGrid<gdst,lwidth> oppsort){
        std::cout<<"wrongfull parent class accsess"<<std::endl;

    }
    void reproduce(std::vector<int> &spotsNeeded, sf::Image lakes_img, sf::Color lake_col){
        std::cout<<"wrongfull parent class accsess"<<std::endl;
    }
    void walk(sf::Image lakes_img, sf::Color lake_col, refGrid<gdst,lwidth> oppsort) {  
        std::cout<<"wrongfull parent class accsess"<<std::endl;
    }
    point<double> toPoint(){
        return point<double>(x,y,id);
    }
    friend std::ostream& operator<<(std::ostream& os, const Creature& m) {
        os << "creature with id: " <<  std::to_string(m.id);
        return os; // enable chaining
    }
    private:
        void draw(sf::RenderTarget& target, sf::RenderStates state) const{

            sf::CircleShape shape(size+age);
            shape.setFillColor(sf::Color::Red);
            shape.setPosition({x-size,y-size});
            target.draw(shape);
            
        }
};

class PreyBasic : public Creature{
    public:
        PreyBasic(){
            
        }
        PreyBasic(double a, double b, bool c, double g){
            age=g;
            x = a;
            y = b;
            isnull = c;
        
        }
        PreyBasic(double a, double b, bool c, double g, int i){
            age=g;
            x = a;
            y = b;
            isnull = c;
            id=i;
        }
    void tick(sf::Image lakes_img, sf::Color lake_col, double *grassArray, std::vector<int> &spotsNeeded, refGrid<gdst,lwidth> oppsort){
        if (!isnull){
            walk(lakes_img, lake_col, oppsort);
            if (hunger <= 0){
                isnull = true;
                return;
            }
            if (age > 1){
                isnull = true;
                return;
            }else if (age > 0.3){
                matetimer -= 0.005 + (0.002*((double)random()/(double)RAND_MAX));
            }

            age+=0.004*(1-hunger)*((double)random()/(double)RAND_MAX);
            hunger -= 0.01;
            eat(oppsort);
            reproduce(spotsNeeded, lakes_img, lake_col);
        
        }
        
    }
    void eat(refGrid<gdst,lwidth> oppsort){
        if (hunger <= 0.3 || iseating){
            iseating = true;
            if (hunger >= 0.9){
                iseating = false;
            }
            int cx = floor(x);
            int cy = floor(y);
            
            if (grassArray[cy*lwidth+cx]-0.1 >= 0){
                grassArray[cy*lwidth+cx]-=0.1;
                hunger += 0.2;
            }else{
                grassArray[cy*lwidth+cx] = 0;
            }
            
        }

    }
    void reproduce(std::vector<int> &spotsNeeded, sf::Image lakes_img, sf::Color lake_col){
        if (matetimer <= 0 && !isnull && !iseating && hunger >= 0.5){
            
            matetimer=1;
            spotsNeeded.push_back(id);

            
        }
    }
    void walk(sf::Image lakes_img, sf::Color lake_col, refGrid<gdst,lwidth> oppsort){
         
  
        if (!isnull){ // check that is exists
            

            //rotate random amount, in radions, "then walk" in that direction
            double temp_rand = (double)random()/(double)RAND_MAX;
            r = (2 * M_PI * temp_rand);

            double movex = cos(r);
            double movey = sin(r);

            x+=movex;
            y+=movey;
            
            if (x < 0){
                x = 0;
            }else 
            if (x >= lwidth){
                x = lwidth-1;
            }
            if (y < 0){
                y = 0;
            }else 
            if (y >= lheight){
                y = lheight-1;
            }
            
            if (lakes_img.getPixel(floor(x),floor(y))==lake_col){
                
                if (lakes_img.getPixel(floor(x),floor(y-movey*0.8))!=lake_col){
                    y-=movey*0.5;
                }else 
                if (lakes_img.getPixel(floor(x-movex*0.5),floor(y))!=lake_col){
                    
                    x-=movex*0.75; 
                }else{
                    x-=movex*0.75;
                    y-=movey*0.75;
                }
            }
        }

        
    }
    private:
        void draw(sf::RenderTarget& target, sf::RenderStates state) const{

            sf::CircleShape shape(size+age);
            shape.setFillColor(sf::Color::Green);
            shape.setPosition({x-size,y-size});
            target.draw(shape);
            
        }

    
};
class PredBasic : public Creature{
    public:
        bool cankill=true;
        int wannakill=-1;
        PredBasic(){
            
        }
        PredBasic(double a, double b, bool c, double g){
            age=g;
            x = a;
            y = b;
            isnull = c;
        
        }
        PredBasic(double a, double b, bool c, double g, int i){
            age=g;
            x = a;
            y = b;
            isnull = c;
            id=i;
        }
    void tick(sf::Image lakes_img, sf::Color lake_col, double *grassArray, std::vector<int> &spotsNeeded, refGrid<gdst,lwidth> oppsort){
        wannakill=-1;
        if (!isnull){
            walk(lakes_img, lake_col, oppsort);
            if (hunger <= 0){
                isnull = true;
                return;
            }
            if (age > 1){
                isnull = true;
                return;
            }else if (age > 0.3){
                matetimer -= 0.006 + (0.002*((double)random()/(double)RAND_MAX));
            }

            age+=0.004*(1-hunger)*((double)random()/(double)RAND_MAX);
            hunger -= 0.01;
            eat(oppsort);
            reproduce(spotsNeeded, lakes_img, lake_col);
        
        }
        
    }
    void eat(refGrid<gdst,lwidth> oppsort){
        if (hunger <= 0.3 || iseating){
            iseating=true;
            point<double> nearopp = oppsort.getNearestAround(toPoint());

            if(nearopp.distance(toPoint()) < 1.2){
                if (hunger+0.2 <= 1){
                    hunger += 0.8;
                }else{
                    hunger = 1;
                }
                wannakill=nearopp.id;
            }
            if (hunger>=0.8){
                iseating=false;
            }
        }

    }
    void reproduce(std::vector<int> &spotsNeeded, sf::Image lakes_img, sf::Color lake_col){
        if (matetimer <= 0 && !isnull && !iseating && hunger >= 0.5){
            
            matetimer=1;
            
            spotsNeeded.push_back(id);

            
        }
    }
    void walk(sf::Image lakes_img, sf::Color lake_col, refGrid<gdst,lwidth> oppsort){  
        if (!isnull){ // check that is exists
    
            //rotate random amount, in radions, "then walk" in that direction
            double temp_rand = (double)random()/(double)RAND_MAX;
            r = (2 * M_PI * temp_rand);

            if (iseating){
                point<double> nearopp = oppsort.getNearestAround(toPoint());
                if (nearopp.id != -1){
                    r = nearopp.direction(toPoint());
                }
            }   

            double movex = cos(r);
            double movey = sin(r);
            
            x+=movex;
            y+=movey;
            
            if (x < 0){
                x = 0;
            }else 
            if (x >= lwidth){
                x = lwidth-1;
            }
            if (y < 0){
                y = 0;
            }else 
            if (y >= lheight){
                y = lheight-1;
            }
            
            if (lakes_img.getPixel(floor(x),floor(y))==lake_col){
                
                if (lakes_img.getPixel(floor(x),floor(y-movey*0.8))!=lake_col){
                    y-=movey*0.5;
                }else 
                if (lakes_img.getPixel(floor(x-movex*0.5),floor(y))!=lake_col){
                    
                    x-=movex*0.75; 
                }else{
                    x-=movex*0.75;
                    y-=movey*0.75;
                }
            }
        }

        
    }
    private:
        void draw(sf::RenderTarget& target, sf::RenderStates state) const{

            sf::CircleShape shape(size+age);
            shape.setFillColor(sf::Color::Red);
            shape.setPosition({x-size,y-size});
            target.draw(shape);
            
        }
 
    
};

// refGrid<PreyBasic 25, lwidth> preysort = refGrid<PreyBasic, 25, lwidth>();
int main(){
    sf::Sprite lakes_spr;
    sf::Sprite grass_spr;


    sf::RenderWindow window(sf::VideoMode(500, 500), "program");
    sf::RenderWindow datawd(sf::VideoMode(500, 500), "data");
    sf::CircleShape shape(100.f);

    sf::Image lakes_img;
    sf::Image grass_img;
    
    
    sf::Color lake_col = sf::Color::Blue;
    
    

    
    
    lakes_img.create(lwidth, lheight, sf::Color{30,5,5});
    grass_img.create(lwidth, lheight, sf::Color{0,0,0,0});
    const siv::PerlinNoise::seed_type seed = 123456u;
    
	const siv::PerlinNoise perlin{ seed };
    
    double lakes_blur[lwidth*lheight];

    

    
    const int GRAPH_MAX = 100;
    int datagraph[GRAPH_MAX];
    
    for (int i = 0; i < GRAPH_MAX; i++){
        datagraph[i]= -1;
    }
    
    for(int y = 0; y < lheight; y++){
        for(int x = 0; x < lwidth; x++){
            
            const double noise = perlin.octave2D_01((x * 0.01), (y * 0.01), 4);
            
            if (noise > 0.45 && noise < 0.6){
                lakes_img.setPixel(x, y, lake_col);
            }else{
                grassArray[y*lwidth+x] = perlin.octave2D_01((x * 0.01), (y * 0.01), 2);
            }
            lakes_blur[y*lwidth+x] = noise;
            
            
        }
    }
    CreatureOrganizer<PreyBasic, CREATURE_MAX> preys = CreatureOrganizer<PreyBasic, CREATURE_MAX>
                        (
                            lakes_img, 
                            lake_col, 
                            CREATURE_START
                        );
    CreatureOrganizer<PredBasic, CREATURE_MAX> preds = CreatureOrganizer<PredBasic, CREATURE_MAX>
                        (
                            lakes_img, 
                            lake_col, 
                            CREATURE_START/2
                        );
    // for (int i = 0; i < CREATURE_MAX; i++){
        
    //     int curx = lwidth * random()/RAND_MAX;
    //     int cury = lheight* random()/RAND_MAX;
    //     while (lakes_img.getPixel(curx, cury) == lake_col){
    //         curx = lwidth * random()/RAND_MAX;
    //         cury = lheight* random()/RAND_MAX;
    //     }

    //     preys[i] = PreyBasic(curx, cury, i >= CREATURE_START, (double)random()/(double)RAND_MAX);

    // }



    sf::Texture lakes_txr;
    

    lakes_txr.loadFromImage(lakes_img);

    lakes_spr.setTexture(lakes_txr);

    int last_octive = 4;

    clock_t deltaTime = 0;
    unsigned int frames = 0;
    double  frameRate = 60;
    double  averageFrameTimeMilliseconds = 33.333;


    
    while (window.isOpen()) {
        clock_t beginFrame = clock();
        
        sf::Event event;
        while (window.pollEvent(event)){
            if (event.type == sf::Event::Closed){
                window.close();

            }
        }
        
        for(int y = 0; y < lheight; y++){
            for(int x = 0; x < lwidth; x++){
                if (lakes_img.getPixel(x,y)==lake_col){
                    grassArray[y*lwidth+x]=0;
                    continue;
                }
               
                int sumcount = 0;
                double near = 0;
                if(x!=0 && lakes_img.getPixel(x-1,y) != lake_col){
                    near += grassArray[y*lwidth+(x-1)];
                    sumcount++;
                }
                if(y!=0 && lakes_img.getPixel(x,y-1) != lake_col){
                    near += grassArray[(y-1)*lwidth+x];
                    sumcount++;
                }
                if(x!=lwidth-1 && lakes_img.getPixel(x+1,y) != lake_col){
                    near += grassArray[y*lwidth+(x+1)];
                    sumcount++;
                }
                if(y!=lheight-1 && lakes_img.getPixel(x,y+1) != lake_col){
                    near += grassArray[(y+1)*lwidth+x];
                    sumcount++;
                }

                if(x!=0 && y!=0 && lakes_img.getPixel(x-1,y-1) != lake_col){
                    near += grassArray[(y-1)*lwidth+(x-1)];
                    sumcount++;
                }
                if(x!=lwidth-1 && y!=0 && lakes_img.getPixel(x+1,y-1) != lake_col){
                    near += grassArray[(y-1)*lwidth+(x+1)];
                    sumcount++;
                }
                if(x!=lwidth-1 && y!=lheight-1 && lakes_img.getPixel(x+1,y+1) != lake_col){
                    near += grassArray[(y+1)*lwidth+(x+1)];
                    sumcount++;
                }
                if(x!=0 && y!=lheight-1 && lakes_img.getPixel(x-1,y+1) != lake_col){
                    near += grassArray[(y+1)*lwidth+(x-1)];
                    sumcount++;
                }
                double nearavg = near/sumcount;
                double incr = 0.0002;
                double cblur = 0.75+(lakes_blur[y*lwidth+x]);
                
                if (nearavg > grassArray[y*lwidth+x]){
                    if (grassArray[y*lwidth+x]+incr*cblur <= 1){
                        grassArray[y*lwidth+x]+=incr*cblur;
                    }else{
                        grassArray[y*lwidth+x] = 1;
                    }
                }
                if (grassArray[y*lwidth+x]-incr*nearavg/cblur >= 0){
                    grassArray[y*lwidth+x]-=incr*nearavg/cblur;
                }else{
                    grassArray[y*lwidth+x] = 0;
                }
                if (grassArray[y*lwidth+x]+incr*(1-(nearavg-grassArray[y*lwidth+x])) <= 1){
                    grassArray[y*lwidth+x]+=incr*(1-(nearavg-grassArray[y*lwidth+x]));
                }else{
                    grassArray[y*lwidth+x] = 1;
                }
                
                
                grass_img.setPixel(x,y,  sf::Color{80, 180, 80, 255 * grassArray[y*lwidth+x]});
            }
        }
        
        window.clear();
        
        sf::Texture grass_txr;
        grass_txr.loadFromImage(grass_img);
        grass_spr.setTexture(grass_txr);
        window.draw(lakes_spr);
        window.draw(grass_spr);

        
        for (int i = 0; i < preytokill.size(); i++){
            preys[preytokill[i]].isnull=true;
        }
        
        preytokill = std::vector<int>();


        preys.tickThreaded(12, preds.gridsort);
        preds.tickThreaded(12, preys.gridsort);
        window.draw(preys);
        window.draw(preds);

        
        if (datawd.isOpen()){
            sf::Event event;
            while (datawd.pollEvent(event)){
                if (event.type == sf::Event::Closed){
                    datawd.close();

                }
            }
            datawd.clear();
            for (int i = 0; i < GRAPH_MAX-1; i++){
                datagraph[i] = datagraph[i+1];
                std::array line = {
                    sf::Vertex{sf::Vector2f(((double)(i-1)/(double)GRAPH_MAX)*500.f, 500.f-500.f*((double)datagraph[i-1]/(double)CREATURE_MAX))},
                    sf::Vertex{sf::Vector2f(((double)i/(double)GRAPH_MAX)*500.f, 500.f-500.f*((double)datagraph[i]/(double)CREATURE_MAX))}
                };
                datawd.draw(line.data(), line.size(), sf::PrimitiveType::Lines);
            }
            datagraph[GRAPH_MAX-1] = preys.pop;
            datawd.display();
            
            
        }
        
        sf::Font font;
        font.loadFromFile("NotoNaskhArabic-Regular.ttf");
        sf::Text fpscount = sf::Text("population: "+std::to_string(preys.pop), font, 15);
        fpscount.setPosition({0,0});
        window.draw(fpscount);

        
        window.display();
        clock_t endFrame = clock();
        deltaTime += endFrame - beginFrame;
        frames++;
        if(clockToMilliseconds(deltaTime)>1000.0){ 
            frameRate = (double)frames*0.5 +  frameRate*0.5; 
            frames = 0;
            deltaTime -= CLOCKS_PER_SEC;
            averageFrameTimeMilliseconds  = 1000.0/(frameRate==0?0.001:frameRate);
        }
    }
    return 0;
}

