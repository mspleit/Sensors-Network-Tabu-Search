#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <fstream>
#include <string>
#include <unordered_map>
#include <map>
#include <algorithm>

#include <time.h>
#include <math.h>

//numerical values of some parameters

const int gridSize = 750; // number of cells of the country, on each side
double pi = 3.14159265358979; //pi
double wind = 15; //wind strength, in km/h
double maxResponseTimeBomb = 4.0; //maximum ellapsed time between attack and its detection
double maxResponseTimePlant = 6.0; //maximum ellapsed time between plant leaks a cloud and its detection
double k_diff = 5; //diffusion term, 
double R_init_attack = 0.3;	//initial radius of the radioactive cloud caused by a terrorist attack (dirty bomb)
							//we model a cloud as a gaussian distriubtion
							//and the chosen radius equals 3*[standard_deviation]
double R_init_plant = 0.9;	//initial radius of the radioactive cloud caused by a plant leaking radioactivity in the air
							//we model a cloud as a gaussian distriubtion
							//and the chosen radius equals 3*[standard_deviation]
double R_init_east = 100; //inital radius of the radioactive cloud coming from the east


						  // costs (in euros) of events
namespace costs {

	double nuclearPlantLeak = 1E+9;
	double easternAccident = 1E+5;
	double dirtyBomb = 1E+6;
	double plantDetected = 0.5 * nuclearPlantLeak;
	double easternDetected = 0.1 * easternAccident;
	double dirtyBombDetected = (2.0 / 3) * dirtyBomb;

	double sensorMaintenance = 5E+3; // maintenance cost per sensor, per year
	double falseAlert = 1E+5; // false alarm cost
	double sensorCheck = 1E+3; // suspicious false alarm check
}

// probabilities (per year)
namespace probs {

	double sensorNotWorking = 1E-6; // a given sensors doesn't work
	double nuclearPlantLeak = 5E-5; // leak from a nuclear plant
	double easternAccident = 1E-4; // accident from the eastern border
	double dirtyBombPerCountry = 1E-5; // dirty bomb in the country
	double dirtyBombPerPixel = dirtyBombPerCountry / (gridSize*gridSize);
	double falseAlarm = 2.0; // false alarm probability (per sensor)
}

// coordinates go from 0 to gridSize-1
struct Point {
	int x, y;
	int dist2(const Point& p) const {
		int dx = p.x - x;
		int dy = p.y - y;
		return dx*dx + dy*dy;
	}
	double dist(const Point& p) const {
		return sqrt(dist2(p));
	}
};

struct Interval;

typedef std::vector<Interval> Intervals;

template<typename T>
T clamp(T x, T minV, T maxV) { return min(maxV, max(minV, x)); }

struct Interval {

	double start, end; // bounds
	int sensor_id; //id of sensor
	Interval(double s = 0, double e = 0, int i = 0) : start(s), end(e), sensor_id(i) {}

	// by default, split in the [0;2*Pi] loop
	Intervals circularSplit(double minV = 0.0, double maxV = 2 * pi) const
	{ // TODO : better PI
	  //always returns two intervals, one of them may be [maxV,maxV]
		double shift = ceil((minV - start) / (maxV - minV))*(maxV - minV);
		double start0 = start + shift,
			end0 = end + shift; // Shift so that "start0" is in [minV,maxV]
		if (end0 <= maxV)
		{
			return{ { start0, end0, sensor_id },{ maxV, maxV, sensor_id } };
		}
		double end2 = minV + (end0 - maxV); // looped end
		if (end2 >= start0)
		{
			// full loop
			return{ { minV, maxV, sensor_id },{ maxV, maxV, sensor_id } };
		}

		return{ { minV, end2, sensor_id },{ start0, maxV, sensor_id } };
	}

	//length of an interval (HACK : assumes start <= end)
	double get_length() const { return end - start; };

	//length an union of intervals. Assume intervals are disjoints and sorted.
	static double get_length(Intervals& inters) {

		if (inters.size() == 0) { return 0; }

		std::sort(inters.begin(), inters.end(),
			[](const Interval& a, const Interval& b) { return a.start <= b.start; });

		double totalLength = 0;

		Interval current = inters[0];

		for (int i = 1; i < inters.size(); i++) {
			const Interval& inter = inters[i];
			if (inter.start > current.end) {
				// no collision -> add current length to total length. Replace current interval
				totalLength += current.get_length();
				current = inter;
			}
			else {
				// collision -> merging i with current interval
				current.end = std::max(current.end, inter.end);
			}
		}

		// add the last interval
		totalLength += current.get_length();

		return totalLength;
	}
};

//define order on intervals
bool operator<(const Interval &a, const Interval &b)
{
	//use lexicographic order on sets
	return ((a.start < b.start) || ((a.start == b.start) && (a.end>b.end)) || ((a.start == b.start) && (a.end == b.end) && (a.sensor_id < b.sensor_id)));
}

// standard deviation (in pixels) of the Eastern accidents
double eastProb( // probability function of the Eastern accidents
	double dist // distance (in pixels) from the eastern border
	) {
	if (dist > 90)
	{
		//assume detection is 0 if distance to border is greater than 90km
		//corresponds to maximum penetration distance in 6 hours with 15km/h wind
		return 0;
	}
	else
	{
		//third order approximation of the probability density
		return 0.3119 - 0.0037*dist + (4E-5f)*(dist*dist) - (3E-7f)*(dist*dist*dist);
	}

	return 0;
}

// standard deviation of a (terrorist) dirty bomb
int bombMaxDist = wind * maxResponseTimeBomb; // distance where bombProp() is approximated by 0
double bombProb( // probability of dectection of a dirty bomb
	double dist // distance (in pixels) from the bomb
	) {
	if (dist > bombMaxDist)
	{
		//no detection is possible if distance is greater than 100km
		//this corresponds to maximum distance reached within 6 hours, wind force 15 km/h
		return 0;
	}
	else
	{
		if (dist <= R_init_attack)
		{

			return 1;
		}
		else
		{
			//always make sure R_init_attack is >0, so you can't divide by zero here
			//analytical value (accurate but takes time to compute)
			return (1 / pi) * acos(((-9 * k_diff / wind) + sqrt((81 * k_diff*k_diff / (wind*wind) + dist*dist - R_init_attack*R_init_attack))) / dist);
			//logarithmic approximation (less accurate but faster)
			//return(0.4323 - 1.2*0.082*log(dist));
		}
	}

	return 0; //should never be reached
}

int plantMaxDist = wind * maxResponseTimePlant;
double plantProb( // probability of dectection of a dirty bomb
	double dist // distance (in pixels) from the bomb
	) {
	if (dist > plantMaxDist || dist <= 5.0)
	{
		//no detection is possible if distance is greater than plantMaxDist
		//this corresponds to maximum distance reached within 6 hours, wind force 15 km/h
		//a quick analysis of the TELERAY network shows that almost no sensor is closer than 5km to a plant
		//EDF's sensors cover 10km around each plant
		return 0;
	}
	else
	{
		if (dist <= R_init_plant)
		{
			//impose no detection if distance is very small, so that sensors located on the plants won't tweak he result
			return 0.0;
		}
		else
		{
			//always make sure R_init_attack is >0, so you can't divide by zero here
			//analytical value (accurate but takes time to compute)
			return (1 / pi) * acos(((-9 * k_diff / wind) + sqrt((81 * k_diff*k_diff / (wind*wind) + dist*dist - R_init_plant*R_init_plant))) / dist);
			//logarithmic approximation (less accurate but faster)
			//return(0.4323 - 1.2*0.082*log(dist));
		}
	}

	return 0; //should never be reached
}

// angular detection interval for a sensor at s, and a pixel at b
Interval computeBombInterval(const Point& b, const Point& s) {

	int dx = s.x - b.x, dy = s.y - b.y;
	double angle = atan2(dy, dx);
	double prob = bombProb(sqrt(dx*dx + dy*dy));
	return{ angle - prob*pi, angle + prob*pi };
};
//compute detection interal for a cloud coming from the east
Interval computeEastInterval(const Point& s)
{
	// eastern border corresponds to x = gridsize-1
	// as detector positions are in [0,gridsize-1]
	int dx = gridSize - 1 - s.x;
	double length = 500.0 * eastProb(dx);
	return{ s.y - length / 2.0, s.y + length / 2.0 };
}
//compute detection interval for one sensor s and one plant p
Interval computePlantInterval(const Point& s, const Point& p)
{
	int dx = p.x - s.x, dy = p.y - s.y;
	double angle = atan2(dy, dx);
	double prob = plantProb(sqrt(dx*dx + dy*dy));
	return{ angle - prob*pi, angle + prob*pi };
}




// coordinates of the nuclear plants
//taking into account the fact that coordinates range from 0 to gridsize-1
std::vector<Point> plants = {

	{ 294, 654 },
	{ 489, 703 },
	{ 693, 698 },
	{ 10, 14 },
	{ 661, 328 },
	{ 472, 432 },
	{ 667, 479 },
	{ 421, 176 },
	{ 360, 261 },
	{ 408, 441 },
	{ 91, 564 },
	{ 632, 33 },
	{ 59, 194 },
	{ 344, 451 },
	{ 264, 650 },
	{ 293, 195 },
	{ 507, 128 },
	{ 280, 460 },
	{ 437, 736 },
	{ 504, 52 }
};



struct SensorNetwork {

private:
	static bool detectionPrecomputed;
	static std::vector<Interval> bombDetection_circular_left;
	static std::vector<Interval> bombDetection_circular_right;
	static std::vector<Interval> plantDetection_circular_left;
	static std::vector<Interval> plantDetection_circular_right;
	static void computeDetectionIntervals()
	{
		//initialize bomb detection intervals
		for (int y = -bombMaxDist; y < bombMaxDist + 1; y++)
		{
			for (int x = -bombMaxDist; x < bombMaxDist + 1; x++)
			{
				Interval inter = computeBombInterval({ x, y }, { 0, 0 });
				Intervals interSplit = inter.circularSplit();
				bombDetection_circular_left.push_back(interSplit[0]);
				bombDetection_circular_right.push_back(interSplit[1]);
			}
		}
		//initialize plant detection intervals
		for (int y = -plantMaxDist; y < plantMaxDist + 1; y++)
		{
			for (int x = -plantMaxDist; x < plantMaxDist + 1; x++)
			{
				Interval inter = computePlantInterval({ x, y }, { 0, 0 });
				Intervals interSplit = inter.circularSplit();
				plantDetection_circular_left.push_back(interSplit[0]);
				plantDetection_circular_right.push_back(interSplit[1]);
			}
		}
	}

public:
	struct Cell {
		double probability = 0; // detected
								// is the probability up to date (with the intervals) ?
		bool updToDate = true;
		std::set<Interval> sensorIntervals;

		void computeProb() {
			double totalLength = 0.0;
			double currentLength = 0.0;
			double currentEnd = 0.0;

			for (const auto& i : sensorIntervals)
			{
				updToDate = true;
				//assume intervals in sensorIntervals are sorted by lexicographic order
				// ie by increasing order start, then end, then id
				if (currentEnd < i.start)
				{
					//no collision : add currentLength, and update current interval
					totalLength += currentLength;
					currentEnd = i.end;
					currentLength = i.end - i.start;
				}
				else
				{
					//collision : add interval i to current interval
					if (i.end > currentEnd)
					{
						currentLength += i.end - currentEnd;
						currentEnd = i.end;

					}
					else
					{
						//current interval already spans interval i
						//therefore : do nothing
					}
				}
				//if current end is 2*pi, no need to go any further
				if (currentEnd == 2 * pi) { break; }
			}

			totalLength += currentLength;

			probability = totalLength / (2 * pi);
		}

		//add an interval to the set
		void add_interval(Interval inter)
		{
			//insertion complies with the order defined on Interval
			//complexity is logarithmic in sensorIntervals size
			if (inter.end > inter.start)
			{
				sensorIntervals.insert(inter);
				updToDate = false;
			}

		}

		//remove a given interval from the set
		void remove_interval(Interval inter)
		{
			if (inter.end > inter.start) //an empty interval is guaranteed not to be in the set
			{
				//complexity is logarithmic in sensorIntervals size
				std::set<Interval>::iterator it = sensorIntervals.find(inter);
				if (it != sensorIntervals.end())
				{
					sensorIntervals.erase(it);
					updToDate = false;
				}
			}
		}

		//replace an interval by another one
		//faster than delete/insert if intervals are close
		void replace_interval(Interval iOld, Interval iNew)
		{
			if (iNew.end == iNew.start)
			{
				//new interval is empty : no need to insert it
				remove_interval(iOld);
				return;
			}
			if (iOld.end == iOld.start)
			{
				//old interval was empty, thus not in the set
				add_interval(iNew);
				return;
			}


			//find old interval
			std::set<Interval>::iterator it = sensorIntervals.find(iOld);
			if (it != sensorIntervals.end())
			{
				sensorIntervals.insert(it, iNew);
				sensorIntervals.erase(it);
				updToDate = false;
			}
			else
			{
				add_interval(iNew);
			}
		}
	};

	struct CellGrid {
		int w, h;
		std::vector<Cell> cells;
		CellGrid(int w, int h) : w(w), h(h), cells(w*h) {};
		Cell& get(int x, int y) {
			//if (x < 0 || x >= w || y < 0 || y >= h) { throw 1; }
			return cells[y*w + x];
		};
	};

	CellGrid cells; // 2D Grid
	std::vector<Cell> plantsCells;	//one cell for each plant
									//plantsCells[i] corresponds to the plant which location is given by plants[i]
	std::map<int, Point> sensors; // all sensors in the network
	std::set<Interval> eastCloud_detectionIntervals;
	double bombTotalBenefits = 0.0; //sum of each cell's individual detection probability, times the associated benefit
	double eastTotalBenefits = 0.0; //total benefits linked to of detecting a radioactive cloud coming from the east
	double maintenanceCosts = 0.0;
	double falseAlarmCosts = 0.0;
	double plantTotalBenefits = 0.0;

	int apprFactor = 1; // approximation factor >= 1, downscales the cell grid (1 for exact computation)

	float epsilon; //maximum distance at which two sensors are still considered "close" to each other
	int epsilonCeil;

	SensorNetwork(int apprFactor = 1, float epsilon = sqrt(2)) :
		apprFactor(apprFactor),
		cells(gridSize / apprFactor, gridSize / apprFactor),
		plantsCells(plants.size()),
		epsilon(epsilon), epsilonCeil(ceil(epsilon))
	{
		// if not done yet, pre-compute the detection intervals
		if (!SensorNetwork::detectionPrecomputed) {
			SensorNetwork::computeDetectionIntervals();
			SensorNetwork::detectionPrecomputed = true;
		}
	}

	// add a new sensor, or change it
	void setSensor(int id, const Point& p) {

		auto oldSensor = sensors.find(id);
		if (oldSensor == sensors.end())
		{
			// if new Id, add it
			sensors[id] = p;

			//update maintenance costs
			maintenanceCosts += costs::sensorMaintenance;

			//update detection for attacks at nearby cells
			Interval iLeft, iRight;

			//add detection Intervals to the newly near cells
			{
				for (int i = std::max(0, p.x - bombMaxDist) / apprFactor; i < std::min(gridSize, p.x + bombMaxDist + 1) / apprFactor; i++)
				{
					for (int j = std::max(0, p.y - bombMaxDist) / apprFactor; j < std::min(gridSize, p.y + bombMaxDist + 1) / apprFactor; j++)
					{
						Cell& cell = cells.get(i, j);

						int x = i*apprFactor + apprFactor / 2;
						int y = j*apprFactor + apprFactor / 2;
						int dx = clamp(x - p.x, -bombMaxDist, bombMaxDist);
						int dy = clamp(y - p.y, -bombMaxDist, bombMaxDist);

						//get the two circularSplit intervals for sensor
						iLeft = bombDetection_circular_left[(bombMaxDist + dx) + (bombMaxDist + dy)*(2 * bombMaxDist + 1)];
						iRight = bombDetection_circular_right[(bombMaxDist + dx) + (bombMaxDist + dy)*(2 * bombMaxDist + 1)];

						//add detection intervals
						//add_interval method garantees that empty intervals are not inserted
						iLeft.sensor_id = id;
						cell.add_interval(iLeft);
						iRight.sensor_id = id;
						cell.add_interval(iRight);

					}
				}
			}

			//update detection for plant leaks
			{
				for (int i = 0; i < plants.size(); i++)
				{
					Cell& cell = plantsCells[i];

					//add detection interval only if plant is in range (otherwise do nothing)
					if (p.dist(plants[i]) <= plantMaxDist)
					{
						int x = plants[i].x;
						int y = plants[i].y;
						int dx = clamp(x - p.x, -plantMaxDist, plantMaxDist);
						int dy = clamp(y - p.y, -plantMaxDist, plantMaxDist);

						//get the two circularSplit intervals for sensor
						iLeft = plantDetection_circular_left[(plantMaxDist + dx) + (plantMaxDist + dy)*(2 * plantMaxDist + 1)];
						iRight = plantDetection_circular_right[(plantMaxDist + dx) + (plantMaxDist + dy)*(2 * plantMaxDist + 1)];

						//add detection intervals
						//add_interval method garantees that empty intervals are not inserted
						iLeft.sensor_id = id;
						cell.add_interval(iLeft);
						iRight.sensor_id = id;
						cell.add_interval(iRight);
					}
				}
			}



		}
		else
		{

			const Point pOld = oldSensor->second;
			if ((pOld.x == p.x) && (pOld.y == p.y)) { return; } //sensor has not moved : no need to do anything

			Interval iLeftOld, iRightOld, iLeftNew, iRightNew;

			//update attack detection probabilities for nearby cells
			{

				//remove Intervals from previously nearby cells
				for (int i = std::max(0, pOld.x - bombMaxDist) / apprFactor; i < std::min(gridSize, pOld.x + bombMaxDist + 1) / apprFactor; i++)
				{
					for (int j = std::max(0, pOld.y - bombMaxDist) / apprFactor; j < std::min(gridSize, pOld.y + bombMaxDist + 1) / apprFactor; j++)
					{
						Cell& cell = cells.get(i, j);

						int x = i*apprFactor + apprFactor / 2;
						int y = j*apprFactor + apprFactor / 2;
						int dx = clamp(x - pOld.x, -bombMaxDist, bombMaxDist);
						int dy = clamp(y - pOld.y, -bombMaxDist, bombMaxDist);

						//get the two circularSplit intervals for sensor
						iLeftOld = bombDetection_circular_left[(bombMaxDist + dx) + (bombMaxDist + dy)*(2 * bombMaxDist + 1)];
						iRightOld = bombDetection_circular_right[(bombMaxDist + dx) + (bombMaxDist + dy)*(2 * bombMaxDist + 1)];

						//erase previous detection intervals
						//remove_interval method guarantees that empty intervals will not be considered, and the up-to-date status is updated correctly
						iLeftOld.sensor_id = id;
						cell.remove_interval(iLeftOld);

						iRightOld.sensor_id = id;
						cell.remove_interval(iRightOld);
					}
				}

				//update Intervals for the newly nearby cells
				for (int i = std::max(0, p.x - bombMaxDist) / apprFactor; i < std::min(gridSize, p.x + bombMaxDist + 1) / apprFactor; i++)
				{
					for (int j = std::max(0, p.y - bombMaxDist) / apprFactor; j < std::min(gridSize, p.y + bombMaxDist + 1) / apprFactor; j++)
					{
						Cell& cell = cells.get(i, j);

						int x = i*apprFactor + apprFactor / 2;
						int y = j*apprFactor + apprFactor / 2;
						int dx = clamp(x - p.x, -bombMaxDist, bombMaxDist);
						int dy = clamp(y - p.y, -bombMaxDist, bombMaxDist);

						//get new detection intervals
						iLeftNew = bombDetection_circular_left[(bombMaxDist + dx) + (bombMaxDist + dy)*(2 * bombMaxDist + 1)];
						iRightNew = bombDetection_circular_right[(bombMaxDist + dx) + (bombMaxDist + dy)*(2 * bombMaxDist + 1)];

						//add detection intervals
						//add_interval method guarantees that empty intervals are not inserted
						//and that the cell's p-to-date status is updated correctly
						iLeftNew.sensor_id = id;
						cell.add_interval(iLeftNew);

						iRightNew.sensor_id = id;
						cell.add_interval(iRightNew);
					}
				}
			}

			//remove old detection interval for eastern cloud
			{
				Interval iEastOld = computeEastInterval(pOld);
				iEastOld.sensor_id = id;
				if (iEastOld.end > iEastOld.start)
				{
					const auto& iter = eastCloud_detectionIntervals.find(iEastOld);
					if (iter != eastCloud_detectionIntervals.end())
					{
						eastCloud_detectionIntervals.erase(iter);
					}
				}
			}
			
			//update detection intervals for plants leaks
			{
				//since loations of all plants are known and fixed, it is not necessary to go through two different cases (ie old and new ranges)
				//both cases can be done one stragiht after the other
				for (int i = 0; i < plants.size(); i++)
				{
					Cell& cell = plantsCells[i];

					int x = plants[i].x;
					int y = plants[i].y;

					//remove detection interval only if plant was in old range (otherwise do nothing)
					if (pOld.dist(plants[i]) <= plantMaxDist)
					{

						int dx = clamp(x - pOld.x, -plantMaxDist, plantMaxDist);
						int dy = clamp(y - pOld.y, -plantMaxDist, plantMaxDist);

						//get the two circularSplit intervals for sensor
						iLeftOld = plantDetection_circular_left[(plantMaxDist + dx) + (plantMaxDist + dy)*(2 * plantMaxDist + 1)];
						iRightOld = plantDetection_circular_right[(plantMaxDist + dx) + (plantMaxDist + dy)*(2 * plantMaxDist + 1)];

						//add detection intervals
						//add_interval method garantees that empty intervals are not inserted
						iLeftOld.sensor_id = id;
						cell.remove_interval(iLeftOld);
						iRightOld.sensor_id = id;
						cell.remove_interval(iRightOld);
					}

					//add new detection intervals only if plant is in new range (otherwise do nothing)
					if (p.dist(plants[i]) <= plantMaxDist)
					{
						int dx = clamp(x - p.x, -plantMaxDist, plantMaxDist);
						int dy = clamp(y - p.y, -plantMaxDist, plantMaxDist);

						//get the two circularSplit intervals for sensor
						iLeftNew = plantDetection_circular_left[(plantMaxDist + dx) + (plantMaxDist + dy)*(2 * plantMaxDist + 1)];
						iRightNew = plantDetection_circular_right[(plantMaxDist + dx) + (plantMaxDist + dy)*(2 * plantMaxDist + 1)];

						//add detection intervals
						//add_interval method garantees that empty intervals are not inserted
						iLeftNew.sensor_id = id;
						cell.add_interval(iLeftNew);
						iRightNew.sensor_id = id;
						cell.add_interval(iRightNew);
					}
				}
			}
			
			sensors[id] = p; // replacing the position
		}

		//TODO : update powerplants detection probabilities and associated costs
		//for reasonnable network size, this can be done explicitely when computing the costs

		//update eastern detection probability
		//add detection interval for eastern cloud, if not empty
		Interval iEastNew = computeEastInterval(p);
		iEastNew.sensor_id = id;
		if (iEastNew.end > iEastNew.start)
		{
			eastCloud_detectionIntervals.insert(iEastNew);
		}

		//TODO (for large sensors networks) : update isolated status
		//otherwise, the computation is done when computing total costs

		return;
	}

	// for each cell, how many sensors are on it ?
	std::vector<int> sensorsCount = std::vector<int>(gridSize*gridSize);

	//compute false alarms costs
	void compute_falseAlarmCosts()
	{
		/*
		//if epsilon is smaller than 20, use the faster version
		if (epsilonCeil < 20)
		{
			// setting the sensorsCount to 0
			for (int i = 0; i < sensorsCount.size(); i++)
			{
				sensorsCount[i] = 0;
			}

			// adding the sensors to their cells

			for (const auto& s : sensors)
			{
				sensorsCount[s.second.y * gridSize + s.second.x]++;
			}

			// counting the number of isolated sensors
			falseAlarmCosts = 0.0;
			for (const auto& s : sensors)
			{
				int x = s.second.x, y = s.second.y;
				bool isolated = true;
				for (int i = std::max(0, x - epsilonCeil); i <= std::min(gridSize - 1, x + epsilonCeil) && isolated; i++)
				{

					for (int j = std::max(0, y - epsilonCeil); j <= std::min(gridSize - 1, y + epsilonCeil) && isolated; j++)
					{
						if ((sensorsCount[y*gridSize + x] > 1) || (sensorsCount[j*gridSize + i] > 0 && ((i - x)*(i - x) + (j - y)*(j - y) <= (epsilon*epsilon)) && ((i != x) || (j != y))))
						{
							isolated = false; break;
						}
					}
				}
				
				falseAlarmCosts += (probs::falseAlarm)*(isolated ? costs::falseAlert : costs::sensorCheck);
			}
		}
		else //epsilon is greater than 20 : use the old fashioned way
		{
			//look at each sensor and see if it is isolated
			falseAlarmCosts = 0.0;
			for (const auto& s : sensors)
			{
				bool isolated = true;
				for (const auto& t : sensors)
				{
					double d = (s.second).dist2(t.second);

					if (((s.first) != (t.first)) && (d <= epsilon*epsilon))
					{
						isolated = false;
						break;
					}
				}
				falseAlarmCosts += (probs::falseAlarm)*(isolated ? costs::falseAlert : costs::sensorCheck);
			}
		}
		*/
		falseAlarmCosts = (probs::falseAlarm)*400*costs::sensorCheck;
	}


	//compute total benefits associated to early detection of eastern clouds
	//value is stored in eastTotalBenefits
	void compute_totalEastBenefits()
	{
		double totalLength = 0.0;
		double currentLength = 0.0;
		double currentEnd = -DBL_MAX;

		for (const auto& i : eastCloud_detectionIntervals)
		{

			//assume intervals in sensorIntervals are sorted by lexicographic order
			// ie by increasing order start, then decreasing end, then increasing id
			if (currentEnd < i.start)
			{
				//no collision : add currentLength, and update current interval
				totalLength += currentLength;
				currentEnd = i.end;
				currentLength = i.end - i.start;
			}
			else
			{
				//collision : add interval i to current interval
				if (i.end > currentEnd)
				{
					currentLength += i.end - currentEnd;
					currentEnd = i.end;

				}
				else
				{
					//current interval already spans interval i
					//therefore : do nothing
				}
			}
		}

		totalLength += currentLength;
		double eastDetectionProb = clamp(totalLength / (gridSize + 500 * eastProb(0.0)), 0.0, 1.0);
		eastTotalBenefits = eastDetectionProb*probs::easternAccident*costs::easternDetected;

		return;
	};

	//total benefits from dirty bombs early detection (all cells are taken into account)
	//value is stored in bombTotalBenefits
	void compute_totalBombBenefits()
	{

		double updateProb = 0; //variation of total detection probability
		for (Cell& cell : cells.cells)
		{
			if (!cell.updToDate)
			{
				cell.computeProb();
				cell.updToDate = true;
			}
			updateProb += cell.probability;
		}

		bombTotalBenefits = (updateProb*apprFactor*apprFactor)*(costs::dirtyBombDetected)*(probs::dirtyBombPerPixel);
		return;
	}

	//compute total benefits associated with detection at power plants
	void compute_totalPlantBenefits()
	{
		int penalty = 0;
		double totalPlantDetectionProb = 0.0;

		//make sure there is at least one sensor at each power plant
		{

			int totalSqrDistance = 0;
			//this vector contains the minimal distance between a plant and the closest sensor (can be zero)
			std::vector<int> dist_min(plants.size(), INT_MAX);

			for (int i = 0; i < plants.size(); i++)
			{
				for (const auto& s : sensors)
				{
					dist_min[i] = std::min(dist_min[i], s.second.dist2(plants[i]));

					if (dist_min[i] == 0) { break; } //minimal distance can't be smaller than 0
				}
			}
			for (int i = 0; i < plants.size(); i++)
			{
				//clamp square distances between 0 and 100, otherwise the total penalty will be greater than INT_MAX (and thus becomes negative)
				//if, for each plant, the closest sensor is more than 100km away, then the total penalty is 2,000,000,000
				//INT_MAX is 2,147,483,647
				totalSqrDistance += clamp(dist_min[i], 0, 100);
			}

			//penalty is 0 if and only if there is a sensor at each plant
			//if so, benefits correspond to early detection of all plants incidents
			penalty = 1000000 * totalSqrDistance;

		}

		//compute detection probabilities for each plant and associated benefits
		{

			for (Cell& cell : plantsCells)
			{
				if (!cell.updToDate)
				{
					cell.computeProb();
					cell.updToDate = true;
				}
				totalPlantDetectionProb += cell.probability;
			}
		}

		plantTotalBenefits = totalPlantDetectionProb*costs::plantDetected*probs::nuclearPlantLeak - penalty;

		return;
	}

	//compute and return total costs
	double getTotalCosts()
	{
		compute_totalBombBenefits();
		compute_totalEastBenefits();
		compute_falseAlarmCosts();
		compute_totalPlantBenefits();

		return maintenanceCosts + falseAlarmCosts - plantTotalBenefits - bombTotalBenefits - eastTotalBenefits;
	}

};

// initialising static variables
bool SensorNetwork::detectionPrecomputed = false;
std::vector<Interval> SensorNetwork::bombDetection_circular_left = std::vector<Interval>();
std::vector<Interval> SensorNetwork::bombDetection_circular_right = std::vector<Interval>();
std::vector<Interval> SensorNetwork::plantDetection_circular_left = std::vector<Interval>();
std::vector<Interval> SensorNetwork::plantDetection_circular_right = std::vector<Interval>();

// testing the network
std::vector<Point> sensors;

// setting random sensors
std::vector<Point> randomSensors(int number = 400) {
	//srand(0);
	std::vector<Point> dst = plants;
	for (int i = plants.size(); i < number; i++) {
		dst.push_back({
			(rand()*(gridSize - 1)) / RAND_MAX,
			(rand()*(gridSize - 1)) / RAND_MAX
		});
	}
	return dst;
}

clock_t start;
void timeStart() { start = clock(); }
void timeEnd() {
	clock_t now = clock();
	std::cout << ((now - start) * 1000) / CLOCKS_PER_SEC << " ms" << std::endl;
	start = now;
}

#include <sstream>
#include <fstream>

using namespace std;
#define	TRUE 1
#define FALSE 0
#define GORANDOM 1

const int num_Sensors = 400;
const int board_size = 750;
// extern std::vector<Point> plants;

class position {
public:
	int x;
	int y;
};

int tabu_list[num_Sensors][board_size][board_size] = { 0 };
int iter;
std::vector<Point> s;
std::vector<std::vector<Point>> neighbors;
std::vector<int> moved;
int best_move;
int step_size;

Point last10[10][num_Sensors];

int memory;
double best_objective;
std::vector<Point> best;

void reinit();
void init();
void get_neighbors(int parent, int step, int startx, int endx, int starty, int endy); //parent is the sensor index
std::vector<Point> get_best();
std::vector<Point> get_best_approx();
bool tabu(int sensor, int x, int y);

SensorNetwork netapprox = SensorNetwork(5,2);
SensorNetwork net = SensorNetwork();

int main(int argc, char **argv)
{
	clock_t begin, end;
	double elapsed_secs, iters_per_sec;
	std::vector<Point> best_neighbor, best_neighbor2;
	bool dynamic_memory = true;
	int i;
	step_size = 1;
	double current_objective = 0;
	double check_objective = 0;
	double iter_best_objective = 0;
	double objective_every_1000 = 0;
	
				 //memory=100;
				 //memory = 4*400; // 4 neighbors per sensor, 400 sensors
				 //memory = 40; // sqrt(4*400)=40
	best_objective = std::numeric_limits<double>::infinity();
	int iters_with_no_change = 0;
	//int max_iter=100000000; //100 million (14s per million = 1400 s ~= 23 mins)
	//int max_iter=10000000; //10 million (14s per million = 140 s ~= 2.5 mins)
	//int max_iter=1.5*60*10; // 1.5 iter/sec*60sec/min * 5 mins ----- 1M (2000 per second = 500s = 8m20s
	int max_iter = 15000; //1 hours
	int search = 10;
	int bigsearch = 200;
	memory = 20;//1500;//250;//board_size; //Initial run

	iter = 1;
	int num_neighbors = 0;
	init(); //initialize positions of Sensors

	double value = net.getTotalCosts(); // = objective(s); // make the initial call to the full objective function

	best_objective = value;
	iter_best_objective = value;

	std::ofstream myfile;
	std::ofstream solsfile;
	std::ofstream currentfile;
	std::string filename = "objective.txt";
	std::string filename2 = "current.txt";
	myfile.open(filename);
	
	std::vector<double> lastobjs, lastobjs2;
	int n = 0;
	begin = clock();
	
	int parent;
	
	while (iter<max_iter)
	{
		end = clock();
		elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		neighbors.clear();
		moved.clear();

		
		srand(time(NULL));
		//parent = rand() % (num_Sensors / 2);
		parent = rand() % (num_Sensors / 2 -  plants.size()) + plants.size();
		//parent++;
		//if (parent >= num_Sensors / 2)
		//	parent = 5 * plants.size();

		int x, y;
		x = s[parent].x;
		y = s[parent].y;

		get_neighbors(parent, 10, 0, 748, 0, 749);
		//get_neighbors(parent, 10, std::max(0, x - bigsearch), std::min(748, x + bigsearch), std::max(0, y - bigsearch), std::min(749, y + bigsearch));
		best_neighbor = get_best_approx();
		int bestx = best_neighbor[parent].x;
		int besty = best_neighbor[parent].y;

		neighbors.clear();
		moved.clear();

		get_neighbors(parent, 1, std::max(0, bestx - search), std::min(748, bestx + search), std::max(0, besty - search), std::min(749, besty + search));
		num_neighbors = neighbors.size();
		best_neighbor = get_best();

		neighbors.clear();
		moved.clear();

		get_neighbors(parent, 1, std::max(0, x - 5), std::min(748, x + 5), std::max(0, y - 5), std::min(749, y + 5));
		best_neighbor2 = get_best();

		neighbors.clear();
		moved.clear();

		neighbors.push_back(best_neighbor);
		neighbors.push_back(best_neighbor2);

		best_neighbor = get_best();

		s[parent].x = best_neighbor[parent].x;
		s[parent].y = best_neighbor[parent].y;
		s[parent + num_Sensors / 2].x = best_neighbor[parent + num_Sensors / 2].x;
		s[parent + num_Sensors / 2].y = best_neighbor[parent + num_Sensors / 2].y;

		net.setSensor(parent, { s[parent].x,s[parent].y });
		net.setSensor(parent+num_Sensors/2, { s[parent + num_Sensors / 2].x,s[parent + num_Sensors / 2].y });
		netapprox.setSensor(parent, { s[parent].x,s[parent].y });
		netapprox.setSensor(parent + num_Sensors / 2, { s[parent + num_Sensors / 2].x,s[parent + num_Sensors / 2].y });
		current_objective = net.getTotalCosts(); //step_objective(s, parent, 1); //best_neighbor changed to s because get_best updates s directly
		/*if (current_objective - iter_best_objective > 1)
		{ 
			for (i = 0; i<num_Sensors; i++)
			{
				//cout << "Sensor " << i << ": " << s[i].x << ", " << s[i].y << endl;
				myfile << s[i].x << ", " << s[i].y << endl;
			}
		}*/
		if (current_objective < iter_best_objective)
		{
			iter_best_objective = current_objective;
		}
		if (current_objective < best_objective)
		{
			for (i = plants.size(); i<num_Sensors; i++)
			{
				best[i].x = s[i].x;
				best[i].y = s[i].y;
			}
			best_objective = current_objective;
		}

		lastobjs2.push_back(iter_best_objective);
		if (lastobjs2.size()>80)
		{
			lastobjs2.erase(lastobjs2.begin());

			if (abs(lastobjs2[lastobjs2.size()] - lastobjs2[0])<0.0001)//If in x iterations the objective hasnt changed much, diversify
			{
				reinit();
				step_size = 1;
				//value = net.getTotalCosts(); 
					//objective(s); //The full objective function needs to be called because we're changing more than 1 neighbor with reinit()
				cout << "Diversifying." << endl;
				iter_best_objective = value;
			}
			lastobjs2.clear();
		}

		if (iter % 1 == 0) {
			cout << setprecision(9) << std::fixed << iter << ", " << current_objective << ", " << best_objective << endl;
			myfile << setprecision(9) << std::fixed << iter << ", " << current_objective << ", " << best_objective << endl;
		}
		if (iter % 50 == 0) {
			end = clock();
			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			iters_per_sec = iter / elapsed_secs;
			cout << setprecision(1) << std::fixed << iters_per_sec << " iterations per second" << endl;

			currentfile.open(filename2);
			for (i = 0; i<num_Sensors; i++)
			{
				//cout << "Sensor " << i << ": " << best[i].x << ", " << best[i].y << endl;
				currentfile << best[i].x << ", " << best[i].y << endl;
			}
			currentfile.close();
		}
		iter++;

	}

	std::cout << "Optimization complete." << endl;

	for (i = 0; i<num_Sensors; i++)
	{
		cout << "Sensor " << i << ": " << best[i].x << ", " << best[i].y << endl;
		myfile << best[i].x << ", " << best[i].y << endl;
	}
	myfile.close();
	std::cin.get(); //wait for keypress
	return 0;
}
void reinit()
{
	int i, j, k;
	Point p;

	if (GORANDOM == 0)
	{
		std::string sensorFileName = "mth6311b.csv";
		//std::string sensorFileName = "500kfh-fix.csv";
		std::ifstream file(sensorFileName.c_str());
		if (!file.is_open()) {
			std::cout << 100000002 << std::endl;
		}
		i = 0;
		while (!file.eof()) {
			std::string line;
			try {
				std::getline(file, line, ',');
				p.x = stoi(line);
				std::getline(file, line);
				p.y = stoi(line);
			}
			catch (...) {
				continue;
			}
			s[i].x = p.x;
			s[i].y = p.y;
			best[i].x = p.x;
			best[i].y = p.y;
			i++;
		}
		file.close();
	}
	//*****************************************

	if (GORANDOM == 1)
	{
		for (i = 0; i < plants.size(); i++)
		{
			p.x = plants[i].x;
			p.y = plants[i].y;
			s.push_back(p);
			best.push_back(p);
		}
		for (i = 0; i < plants.size(); i++) {
			p.x = plants[i].x - 7;
			p.y = plants[i].y;
			s.push_back(p);
			best.push_back(p);
		}
		for (i = 0; i < plants.size(); i++) {
			p.x = plants[i].x;
			p.y = plants[i].y + 6;
			s.push_back(p);
			best.push_back(p);
		}
		for (i = 0; i < plants.size(); i++) {
			p.x = plants[i].x + 6;
			p.y = plants[i].y;
			s.push_back(p);
			best.push_back(p);
		}
		for (i = 0; i < plants.size(); i++) {
			p.x = plants[i].x;
			p.y = plants[i].y - 6;
			s.push_back(p);
			best.push_back(p);
		}

		srand(time(NULL));
		for (i = plants.size() * 5; i < (num_Sensors) / 2; i++) //Init the first half of non-plant sensors
		{
			p.x = rand() % (board_size-1);
			p.y = rand() % board_size;
			s[i].x = p.x;
			s[i].y = p.y;
		}
		for (i = (num_Sensors) / 2; i < num_Sensors; i++) //Init twins of the first half of sensors
		{
			p.x = s[i - (int)(num_Sensors) / 2].x + 1; //Offset the twin to the right
			//if (p.x == board_size) //Put the twin on the left if we're at the right edge of the board
			//	p.x = board_size - 1;
			p.y = s[i - (int)(num_Sensors) / 2].y;
			s[i].x = p.x;
			s[i].y = p.y;
		}
		for (i = 0; i < num_Sensors; i++)
		{
			net.setSensor(i, { s[i].x, s[i].y });
			netapprox.setSensor(i, { s[i].x, s[i].y });
		}
	}
}

void init()
{
	int i, j, k;
	Point p;

	if (GORANDOM == 1)
	{
		//*****************************************

		for (i = 0; i < plants.size(); i++)
		{
			p.x = plants[i].x;
			p.y = plants[i].y;
			s.push_back(p);
			best.push_back(p);
		}
		for (i = 0; i < plants.size(); i++){
			p.x = plants[i].x - 7;
			p.y = plants[i].y;
			s.push_back(p);
			best.push_back(p);
		}
		for (i = 0; i < plants.size(); i++){
			p.x = plants[i].x;
			p.y = plants[i].y + 6;
			s.push_back(p);
			best.push_back(p);
		}
		for (i = 0; i < plants.size(); i++){
			p.x = plants[i].x + 6;
			p.y = plants[i].y;
			s.push_back(p);
			best.push_back(p);
		}
		for (i = 0; i < plants.size(); i++){
			p.x = plants[i].x;
			p.y = plants[i].y - 6;
			s.push_back(p);
			best.push_back(p);
		}
		
		/*for(i=plants.size(); i<6; i++)
		{  p.x = 749;
		p.y =rand() % board_size +1;
		}*/
		srand(time(NULL));
		for (i = plants.size()*5; i<(num_Sensors) / 2; i++) //Init the first half of non-plant sensors
		{
			p.x = rand() % (board_size-1);
			p.y = rand() % board_size;
			s.push_back(p);
			best.push_back(p);
		}
		for (i = (num_Sensors) / 2; i< num_Sensors; i++) //Init twins of the first half of sensors
		{
			p.x = s[i - (int)(num_Sensors) / 2].x + 1; //Offset the twin to the right
			//if (p.x == board_size) //Put the twin on the left if we're at the right edge of the board
			//	p.x = board_size - 1;
			p.y = s[i - (int)(num_Sensors) / 2].y;
			s.push_back(p);
			best.push_back(p);
		}
	
	}
	//*****************************************
	// Reading the sensors file
	if (GORANDOM == 0)
	{
		std::string sensorFileName = "mth6311b.csv";
		//std::string sensorFileName = "500kfh-fix.csv";
		std::ifstream file(sensorFileName.c_str());
		if (!file.is_open()) {
			std::cout << 100000002 << std::endl;
		}

		while (!file.eof()) {
			std::string line;
			try {
				std::getline(file, line, ',');
				p.x = stoi(line);
				std::getline(file, line);
				p.y = stoi(line);
			}
			catch (...) {
				continue;
			}
			s.push_back(p);
			best.push_back(p);
		}
		file.close();
	}
	for (i = 0; i<num_Sensors; i++)
	{
		net.setSensor(i, { s[i].x, s[i].y });
		netapprox.setSensor(i, { s[i].x, s[i].y });
	}
}

void get_neighbors(int parent, int step, int startx, int endx, int starty, int endy) //parent is the sensor index
{
	int i, j, k;
	vector<Point> p;
	Point q;
	double obj;
	int x, y;

	//Create a new base neighbor solution. We will only update the moving sensor to create children.
	for (k = 0; k<num_Sensors; k++)
	{
		q.x = s[k].x;
		q.y = s[k].y;
		p.push_back(q);
	}

	srand(time(NULL));
	for (i = startx; i<endx; i += step)
	{
		for (j = starty; j<endy; j += step)
		{
			if (step>1)
			{
				x = i + (rand() % (2 * step)) - step + step / 2;
				y = j + (rand() % (2 * step)) - step + step / 2;
			}
			else
			{
				x = i;
				y = j;
			}
			if (x<0) x = 0;
			if (x>748) x = 748;
			if (y<0) y = 0;
			if (y>749) y = 749;

			if (!tabu(parent, x, y) && !(x == s[parent].x && y == s[parent].y))
			{ //if the step isn't Tabu, make it a viable neighbor
				p[parent].x = x;
				p[parent].y = y;

				p[parent + num_Sensors / 2].x = x + 1;
				p[parent + num_Sensors / 2].y = y;

				neighbors.push_back(p);
				moved.push_back(parent);
			}
			else
			{ //if the step IS Tabu, check if it's the potentially the optimum, if it is, update the optimum
				int oldx = s[parent].x;
				int oldy = s[parent].y;
				int oldx2 = s[parent + num_Sensors / 2].x;
				int oldy2 = s[parent + num_Sensors / 2].y;

				s[parent].x = x;
				s[parent].y = y;

				s[parent + num_Sensors / 2].x = x + 1; //check for going off the board!!!!
				s[parent + num_Sensors / 2].y = y;
				
				net.setSensor(parent, { s[parent].x, s[parent].y });
				net.setSensor(parent + num_Sensors / 2, { s[parent + num_Sensors / 2].x, s[parent + num_Sensors / 2].y });
				obj = net.getTotalCosts();

				//obj = step_objective(s, parent, 0);
				//obj=objective(s);
				if (obj<best_objective) {
					best[parent].x = s[parent].x;
					best[parent].y = s[parent].y;
					best[parent + num_Sensors / 2].x = s[parent + num_Sensors / 2].x;
					best[parent + num_Sensors / 2].y = s[parent + num_Sensors / 2].y;
					best_objective = obj;
				}
				s[parent].x = oldx;
				s[parent].y = oldy;
				s[parent + num_Sensors / 2].x = oldx2;
				s[parent + num_Sensors / 2].y = oldy2;
				//net.setSensor(parent, { s[parent].x, s[parent].y });
				//net.setSensor(parent + num_Sensors / 2, { s[parent + num_Sensors / 2].x, s[parent + num_Sensors / 2].y });
			}
		}
	}
}


vector<Point> get_best()
{
	int max;
	int i, k;
	double value;
	double best_neighbor_value = std::numeric_limits<double>::infinity();;
	int best_neighbor_index = 0;
	vector<Point> s_prime;
	Point p;

	for (k = 0; k<num_Sensors; k++)
	{
		p.x = s[k].x;
		p.y = s[k].y;
		s_prime.push_back(p);
	}
	max = neighbors.size();

	for (i = 0; i<max; i++)
	{
		int a = neighbors[i].size();
		
		net.setSensor(moved[i], { neighbors[i][moved[i]].x, neighbors[i][moved[i]].y });
		net.setSensor(moved[i] + num_Sensors / 2, { neighbors[i][moved[i] + num_Sensors / 2].x, neighbors[i][moved[i] + num_Sensors / 2].y });
		value = net.getTotalCosts();
		//net.setSensor(moved[i], { s[moved[i]].x, s[moved[i]].y });
		//net.setSensor(moved[i] + num_Sensors / 2, { s[moved[i] + num_Sensors / 2].x, s[moved[i] + num_Sensors / 2].y });

		//value = step_objective(neighbors[i], moved[i], 0);
		//value=objective(neighbors[i]);
		if (value < best_neighbor_value)
		{
			best_neighbor_value = value;
			best_neighbor_index = i;
			best_move = moved[i];
		}
	}
	tabu_list[moved[best_neighbor_index]][s[moved[best_neighbor_index]].x][s[moved[best_neighbor_index]].y] = iter + memory;
	//iterate through neighbors
	//for each, calculate objective
	//keep track of best objective
	//whichever Sensor is the one moving, make the old x,y of that Sensor tabu
	//at the end, return the best neighbor
	/*
	for (i=0; i<num_Sensors; i++)
	{
	p.x = neighbors[best_neighbor_index][i].x;
	p.y = neighbors[best_neighbor_index][i].y;
	s_prime.push_back(p);

	s[i].x = neighbors[best_neighbor_index][i].x;
	s[i].y = neighbors[best_neighbor_index][i].y;
	}*/
	//return s_prime;
	/*
	s[best_move].x=neighbors[best_neighbor_index][best_move].x;
	s[best_move].y=neighbors[best_neighbor_index][best_move].y;
	s[best_move+num_Sensors/2].x=neighbors[best_neighbor_index][best_move].x;
	s[best_move+num_Sensors/2].y=neighbors[best_neighbor_index][best_move].y;
	*/
	s_prime[best_move].x = neighbors[best_neighbor_index][best_move].x;
	s_prime[best_move].y = neighbors[best_neighbor_index][best_move].y;
	s_prime[best_move + num_Sensors / 2].x = neighbors[best_neighbor_index][best_move + num_Sensors / 2].x;
	s_prime[best_move + num_Sensors / 2].y = neighbors[best_neighbor_index][best_move + num_Sensors / 2].y;

	return s_prime;
}
vector<Point> get_best_approx()
{
	int max;
	int i, k;
	double value;
	double best_neighbor_value = std::numeric_limits<double>::infinity();;
	int best_neighbor_index = 0;
	vector<Point> s_prime;
	Point p;

	for (k = 0; k<num_Sensors; k++)
	{
		p.x = s[k].x;
		p.y = s[k].y;
		s_prime.push_back(p);
	}
	max = neighbors.size();

	for (i = 0; i<max; i++)
	{
		int a = neighbors[i].size();

		netapprox.setSensor(moved[i], { neighbors[i][moved[i]].x, neighbors[i][moved[i]].y });
		netapprox.setSensor(moved[i] + num_Sensors / 2, { neighbors[i][moved[i] + num_Sensors / 2].x, neighbors[i][moved[i] + num_Sensors / 2].y });
		value = netapprox.getTotalCosts();
		//net.setSensor(moved[i], { s[moved[i]].x, s[moved[i]].y });
		//net.setSensor(moved[i] + num_Sensors / 2, { s[moved[i] + num_Sensors / 2].x, s[moved[i] + num_Sensors / 2].y });

		//value = step_objective(neighbors[i], moved[i], 0);
		//value=objective(neighbors[i]);
		if (value < best_neighbor_value)
		{
			best_neighbor_value = value;
			best_neighbor_index = i;
			best_move = moved[i];
		}
	}
	tabu_list[moved[best_neighbor_index]][s[moved[best_neighbor_index]].x][s[moved[best_neighbor_index]].y] = iter + memory;
	//iterate through neighbors
	//for each, calculate objective
	//keep track of best objective
	//whichever Sensor is the one moving, make the old x,y of that Sensor tabu
	//at the end, return the best neighbor
	/*
	for (i=0; i<num_Sensors; i++)
	{
	p.x = neighbors[best_neighbor_index][i].x;
	p.y = neighbors[best_neighbor_index][i].y;
	s_prime.push_back(p);

	s[i].x = neighbors[best_neighbor_index][i].x;
	s[i].y = neighbors[best_neighbor_index][i].y;
	}*/
	//return s_prime;
	/*
	s[best_move].x=neighbors[best_neighbor_index][best_move].x;
	s[best_move].y=neighbors[best_neighbor_index][best_move].y;
	s[best_move+num_Sensors/2].x=neighbors[best_neighbor_index][best_move].x;
	s[best_move+num_Sensors/2].y=neighbors[best_neighbor_index][best_move].y;
	*/
	s_prime[best_move].x = neighbors[best_neighbor_index][best_move].x;
	s_prime[best_move].y = neighbors[best_neighbor_index][best_move].y;
	s_prime[best_move + num_Sensors / 2].x = neighbors[best_neighbor_index][best_move + num_Sensors / 2].x;
	s_prime[best_move + num_Sensors / 2].y = neighbors[best_neighbor_index][best_move + num_Sensors / 2].y;

	return s_prime;
}


bool tabu(int sensor, int x, int y)
{
	int i;
	bool b;
	i = tabu_list[sensor][x][y];
	if (tabu_list[sensor][x][y] < iter)
		b = false;
	else
		b = true;
	return b;
}

