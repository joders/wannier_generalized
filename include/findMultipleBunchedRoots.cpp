#ifndef _FIND_MULTIPLE_BUNCHED_ROOTS_CPP_
#define _FIND_MULTIPLE_BUNCHED_ROOTS_CPP_

#include<set>
#include<cfloat>
#include<cmath>
#include<iostream>
#include<iomanip>
#include"jHelper.cpp"
#include"utils.cpp"

using namespace std;

class findMultipleBunchedRoots {
    public:
        mt intervalStart;             // best results, if this is set slightly below ground state energy
        mt intervalEnd;
        mt closestRootSpacing;
        mt (*f)(mt);

        mt gapFactorToNewBunch;

        mt iterationsBeforeStepIncrementation;
        mt incrementationFactor;
        mt stepsBetweenRootsInBunch;  // amount of steps between roots when roots were to be equidistantly spaced

        mt rootValueAccuracy;

        mt rangeDivision;

        set<mt> allRoots;

        findMultipleBunchedRoots(): intervalStart(0.), intervalEnd(100.), closestRootSpacing(.001), f(NULL),
                                    gapFactorToNewBunch(2.), 
                                    iterationsBeforeStepIncrementation(32), incrementationFactor(2.), stepsBetweenRootsInBunch(16), 
                                    rootValueAccuracy(1E-12), 
                                    rangeDivision(0.) {}

        void searchDownUp() {
            searchDownFromInsideBunchesDownwardsThenUpOnEveryMidBunchGapEdge(intervalEnd,intervalStart);
        }
        void searchUpDown() {
            searchUpFromInsideBunchesUpwardsThenDownOnEveryMidBunchGapEdge(intervalStart,intervalEnd);
        }

        void searchZigZag() {

/*
             search works like this:
        
           stopValue         /\
                            /  \
           rangeDivision   /    \
                                 \  /
           startValue             \/
        
                           t_0      t_end           
*/

            //rangeDivision=(4*intervalStart+4*intervalEnd)/8.;
            cout ex(rangeDivision) << endl;

            searchUpFromInsideBunchesUpwardsThenDownOnEveryMidBunchGapEdge(rangeDivision, intervalEnd);
            searchDownFromInsideBunchesDownwardsThenUpOnEveryMidBunchGapEdge(rangeDivision, intervalStart);
        }

        struct DEBUG {
            static const bool PROPAGATE_CALLS, SEARCH_CALLS, FOUND_ROOT; };

    private:
        void searchUpFromInsideBunchesUpwardsThenDownOnEveryMidBunchGapEdge(mt const lowerBorder, mt const upperBorder) {
            foundRootSequence lastRoots(gapFactorToNewBunch);
            searchFromWithinBunches(lowerBorder, +closestRootSpacing, upperBorder, lastRoots);
            lastRoots.reset();
            for(set<mt>::reverse_iterator it =set<mt>::reverse_iterator(allRoots.upper_bound(upperBorder));
                                          it!=set<mt>::reverse_iterator(allRoots.upper_bound(lowerBorder));) {
                lastRoots.pushRoot(*it);
                foundRootSequence gapTester(lastRoots, gapFactorToNewBunch);
                it++;
                if(it!=allRoots.rend()) {
                    gapTester.pushRoot(*it);
                    if((gapTester.thisRootStatus()==gapTester.isInNewBunch) ||
                       (gapTester.thisRootStatus()==gapTester.secondAdded)    ) {

                        searchFromWithinBunches(lastRoots.thisRoot-closestRootSpacing, -closestRootSpacing, lowerBorder/**it+closestRootSpacing*/, lastRoots);
                        it=set<mt>::reverse_iterator(++allRoots.lower_bound(gapTester.thisRoot)); // reverse_iterator constructor shifts one position
                    }
                } else {
                    searchFromWithinBunches(lastRoots.thisRoot-closestRootSpacing, -closestRootSpacing, lowerBorder, lastRoots);
                    break;
                }
            }
        }

        void searchDownFromInsideBunchesDownwardsThenUpOnEveryMidBunchGapEdge(mt const upperBorder, mt const lowerBorder) {
            foundRootSequence lastRoots(gapFactorToNewBunch);
            searchFromWithinBunches(upperBorder, -closestRootSpacing, lowerBorder, lastRoots);
            lastRoots.reset();
            for(set<mt>::iterator it =allRoots.lower_bound(lowerBorder);
                                  it!=allRoots.lower_bound(upperBorder);) {
                lastRoots.pushRoot(*it);
                foundRootSequence gapTester(lastRoots, gapFactorToNewBunch);
                it++;
                if(it!=allRoots.end()) {
                    gapTester.pushRoot(*it);
                    if((gapTester.thisRootStatus()==gapTester.isInNewBunch) ||
                       (gapTester.thisRootStatus()==gapTester.secondAdded)    ) {
                        searchFromWithinBunches(lastRoots.thisRoot+closestRootSpacing, +closestRootSpacing, upperBorder/**it-closestRootSpacing*/, lastRoots);
                        it=allRoots.lower_bound(gapTester.thisRoot);
                    }
                } else {
                    searchFromWithinBunches(lastRoots.thisRoot+closestRootSpacing, +closestRootSpacing, upperBorder, lastRoots);
                    break;
                }
            }
        }

    public:
        struct foundRootSequence { // stores the last three Es that are pushed into it
            mt gapFactorToNewBunch;
            foundRootSequence(mt const gapFactorToNewBunch=2.) {
                reset(); 
                this->gapFactorToNewBunch=gapFactorToNewBunch;
            }
            foundRootSequence(foundRootSequence const &toCopy, mt const gapFactorToNewBunch=2.) {
                thisRoot=toCopy.thisRoot;
                previousRoot=toCopy.previousRoot;
                prepreviousRoot=toCopy.prepreviousRoot;
                this->gapFactorToNewBunch=gapFactorToNewBunch;
            }
            mt thisRoot, previousRoot, prepreviousRoot;
            void pushRoot(mt const arg) {
                prepreviousRoot=previousRoot;
                previousRoot=thisRoot;
                thisRoot=arg;
            }
            void reset() {
                thisRoot=DBL_MAX;
                previousRoot=DBL_MAX;
                prepreviousRoot=DBL_MAX;
            }
            enum rootStatus {isInNewBunch, isInOldBunch, noneAdded, firstAdded, secondAdded};
            rootStatus thisRootStatus() const {
                if(thisRoot==DBL_MAX)
                    return noneAdded;
                else if(previousRoot==DBL_MAX)
                    return firstAdded;
                else if(prepreviousRoot==DBL_MAX)
                    return secondAdded;
                else if(fabs(thisRoot-previousRoot)>gapFactorToNewBunch*fabs(previousRoot-prepreviousRoot))
                    return isInNewBunch;
                else
                    return isInOldBunch;
            }
            mt thisGap() const {
                if( (thisRoot==DBL_MAX) ||
                    (previousRoot==DBL_MAX) )
                    return DBL_MAX;

                return(thisRoot-previousRoot);
            }
        };

    private:
        mt step;
        void searchFromWithinBunches(mt const startFrom, mt const initStep, mt const endWith, foundRootSequence &lastRoots) {
            if(DEBUG::SEARCH_CALLS) cout << "searchFromWithinBunches(" << startFrom << ", " << initStep << ", " << endWith << ", lastRoots)" << endl;
            bool foundNewRoot;
            nonZeroOrderedPair formerBoundary(f,rootValueAccuracy);
            nonZeroOrderedPair laterBoundary(f,rootValueAccuracy);

            step=initStep;
            laterBoundary.setX(startFrom);
            do {
                int itCount=1;
                do {
                    formerBoundary=laterBoundary;
                    if(isAfterEnd(formerBoundary.getX(), initStep, endWith))
                        break;
                    laterBoundary.setX(laterBoundary.getX()+step);
                    if(haveDifferentSign(formerBoundary.getY(),laterBoundary.getY()))
                        break;
                    if(DEBUG::PROPAGATE_CALLS) cout << "P"; cout.flush();
                    if(itCount>=iterationsBeforeStepIncrementation) {
                        step*=incrementationFactor;
                        itCount=0;
                    }
                    itCount++;
                } while(true);
                if(isAfterEnd(formerBoundary.getX(), initStep, endWith))
                    break;

                mt foundRoot;
                if(formerBoundary.getX()<laterBoundary.getX())
                    foundRoot=Ridder_find_val(f,0.,formerBoundary.getX(),laterBoundary.getX(),rootValueAccuracy);
                else
                    foundRoot=Ridder_find_val(f,0.,laterBoundary.getX(),formerBoundary.getX(),rootValueAccuracy);

                foundNewRoot=addNewRoot(foundRoot,lastRoots);

                adjustStepSize(initStep,lastRoots);
                laterBoundary.setX(foundRoot+initStep);

                if(DEBUG::FOUND_ROOT) {
                    if(lastRoots.thisRootStatus()==lastRoots.isInNewBunch)
                        cout << "GAP: " << fabs(lastRoots.previousRoot-lastRoots.thisRoot) << endl;
                    cout ex(foundRoot) << endl;
                }

            } while(foundNewRoot==true);
        };

        class nonZeroOrderedPair {
            mt x,y;
            mt (*f)(mt);
            mt advance;
            public:
                nonZeroOrderedPair(mt (*f_init)(mt), mt advance=1E-12) { f=f_init; this->advance=advance; }
                void setX(mt arg) {
                    y=0;
                    while(true) {
                        y=f(arg);
                        if(y!=0)
                            break;
                        arg=arg*(1.+advance)+advance;
                    }
                    x=arg;
                }
                mt getX() { return x; }
                mt getY() { return y; }
        };

        bool addNewRoot(mt const root, foundRootSequence &lastRoots) {
            set<mt>::iterator it=allRoots.lower_bound(root);
            bool rootPresent=false;
            if(it!=allRoots.end())
                if(fabs(*it-root)<=2.*rootValueAccuracy)
                    rootPresent=true;
            if(it!=allRoots.begin()) {
                it--;
                if(fabs(*it-root)<=2.*rootValueAccuracy)
                    rootPresent=true;
            }

            if(rootPresent==true)
                return false;
            else {
                allRoots.insert(root);
                lastRoots.pushRoot(root);

                return true;
            }
        }

        void adjustStepSize(mt const initStep, foundRootSequence const &lastRoots) {
            if((lastRoots.thisRootStatus()==lastRoots.isInNewBunch) ||
               (lastRoots.thisRootStatus()==lastRoots.firstAdded)  ||
               (lastRoots.thisRootStatus()==lastRoots.secondAdded)    )
                step=initStep;
            else if(lastRoots.thisRootStatus()==lastRoots.isInOldBunch)
                step=lastRoots.thisGap()/stepsBetweenRootsInBunch;
        }

        bool isBeforeEnd(mt const laterBoundary, mt const initStep, mt const endWith) const {
            if(initStep>0)
                return(laterBoundary<endWith);
            else
                return(laterBoundary>endWith);
        }
        bool isAfterEnd(mt const laterBoundary, mt const initStep, mt const endWith) const {
            return(!isBeforeEnd(laterBoundary, initStep, endWith)); }

        bool isPos(mt const arg) const {
            if(arg>=0)
                return(true);
            else
                return(false);
        }
        bool isNeg(mt const arg) const {
            return(!isPos(arg));
        }
        bool haveSameSign(mt const arg1, mt const arg2) const {
            return((isPos(arg1)&&isPos(arg2)) || (isNeg(arg1)&&isNeg(arg2)));
        }
        bool haveDifferentSign(mt const arg1, mt const arg2) const {
            return(!haveSameSign(arg1,arg2));
        }

};

#endif // _FIND_MULTIPLE_BUNCHED_ROOTS_CPP_
