#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cstdlib>

using namespace std;

double NEG_INF = -1e300;

double safe_log(double p) {
    if (p <= 0.0) return NEG_INF;
    return log(p);
}

void trim(string &s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a==string::npos) { s=""; return; }
    size_t b = s.find_last_not_of(" \t\r\n");
    s = s.substr(a, b-a+1);
}

bool starts_with(const string &s, const string &pref) {
    if (s.size() < pref.size()) return false;
    return s.compare(0, pref.size(), pref) == 0;
}

int find_index(string *arr, int n, const string &key) {
    for (int i=0;i<n;++i) if (arr[i] == key) return i;
    return -1;
}

string *split_to_array(const string &s, int &out_n) {
    istringstream iss(s);
    string tok;
    out_n = 0;
    while (iss >> tok) ++out_n;
    if (out_n == 0) return NULL;
    string *arr = new string[out_n];
    istringstream iss2(s);
    int i=0;
    while (iss2 >> tok) { arr[i++] = tok; }
    return arr;
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string model_file;
    string obs_file;
    if (argc == 3) {
        model_file = argv[1];
        obs_file = argv[2];
    } else if (argc == 2) {
        model_file = argv[1];
        obs_file = "";
    } else {
        cerr << "Usage: " << argv[0] << " model.txt [observations.txt]\n";
        return 1;
    }

    ifstream fin(model_file);
    if (!fin) { cerr << "Cannot open model file\n"; return 1; }

    enum Section { NONE, START, TRANS, EMIT };
    Section sec = NONE;

    string *state_names = NULL;
    int S = 0;
    string *symbol_names = NULL;
    int M = 0;

    double *start_probs = NULL;
    double **trans_probs = NULL;
    double **emit_probs = NULL;

    string line;
    while (getline(fin, line)) {
        trim(line);
        if (line.empty()) continue;
        if (starts_with(line, "states:")) {
            string rest = line.substr(7);
            trim(rest);
            if (state_names) delete [] state_names;
            state_names = split_to_array(rest, S);
            if (S <= 0) { cerr<<"states list empty\n"; return 1; }
            start_probs = new double[S];
            for (int i=0;i<S;++i) start_probs[i]=0.0;
            trans_probs = new double*[S];
            for (int i=0;i<S;++i) {
                trans_probs[i] = new double[S];
                for (int j=0;j<S;++j) trans_probs[i][j]=0.0;
            }
            continue;
        }
        if (starts_with(line, "symbols:")) {
            string rest = line.substr(8);
            trim(rest);
            if (symbol_names) delete [] symbol_names;
            symbol_names = split_to_array(rest, M);
            if (M <= 0) { cerr<<"symbols list empty\n"; return 1; }
            emit_probs = new double*[S];
            for (int i=0;i<S;++i) {
                emit_probs[i] = new double[M];
                for (int k=0;k<M;++k) emit_probs[i][k]=0.0;
            }
            continue;
        }
        if (starts_with(line, "start:")) { sec = START; continue; }
        if (starts_with(line, "trans:")) { sec = TRANS; continue; }
        if (starts_with(line, "emit:")) { sec = EMIT; continue; }

        if (sec == START) {
            istringstream iss(line);
            string s; double p;
            if (!(iss >> s >> p)) { cerr<<"Bad start line\n"; return 1; }
            int idx = find_index(state_names, S, s);
            if (idx == -1) { cerr<<"Unknown state in start\n"; return 1; }
            start_probs[idx] = p;
        } else if (sec == TRANS) {
            istringstream iss(line);
            string a,b; double p;
            if (!(iss >> a >> b >> p)) { cerr<<"Bad trans line\n"; return 1; }
            int ia = find_index(state_names, S, a);
            int ib = find_index(state_names, S, b);
            if (ia==-1 || ib==-1) { cerr<<"Unknown state in trans\n"; return 1; }
            trans_probs[ia][ib] = p;
        } else if (sec == EMIT) {
            istringstream iss(line);
            string st, sym; double p;
            if (!(iss >> st >> sym >> p)) { cerr<<"Bad emit line\n"; return 1; }
            int is = find_index(state_names, S, st);
            int ik = find_index(symbol_names, M, sym);
            if (is==-1 || ik==-1) { cerr<<"Unknown state/symbol in emit\n"; return 1; }
            emit_probs[is][ik] = p;
        } else {
            // ignore lines outside sections
        }
    }
    fin.close();

    if (!state_names || !symbol_names) { cerr<<"Model must declare states and symbols (first)\n"; return 1; }

    string *obs = NULL;
    int T = 0;
    if (obs_file != "") {
        ifstream ofn(obs_file);
        if (!ofn) { cerr<<"Cannot open observations file\n"; return 1; }
        string sline;
        if (!getline(ofn, sline)) { cerr<<"Observations file empty\n"; return 1; }
        obs = split_to_array(sline, T);
        ofn.close();
    } else {
        string sline;
        if (!getline(cin, sline)) { cerr<<"No observations on stdin\n"; return 1; }
        obs = split_to_array(sline, T);
    }
    if (T <= 0) { cerr<<"No observations\n"; return 1; }

    int *obs_idx = new int[T];
    for (int t=0;t<T;++t) {
        int id = find_index(symbol_names, M, obs[t]);
        if (id == -1) { cerr<<"Unknown observation symbol\n"; return 1; }
        obs_idx[t] = id;
    }

    double *log_start = new double[S];
    double **log_trans = new double*[S];
    double **log_emit = new double*[S];
    for (int i=0;i<S;++i) {
        log_start[i] = safe_log(start_probs[i]);
        log_trans[i] = new double[S];
        log_emit[i] = new double[M];
        for (int j=0;j<S;++j) log_trans[i][j] = safe_log(trans_probs[i][j]);
        for (int k=0;k<M;++k) log_emit[i][k] = safe_log(emit_probs[i][k]);
    }

    int **backptr = new int*[T];
    for (int t=0;t<T;++t) {
        backptr[t] = new int[S];
        for (int s=0;s<S;++s) backptr[t][s] = -1;
    }

    double *dp_prev = new double[S];
    double *dp_curr = new double[S];

    for (int s=0;s<S;++s) dp_prev[s] = log_start[s] + log_emit[s][obs_idx[0]];

    for (int t=1;t<T;++t) {
        for (int cur=0; cur<S; ++cur) {
            double emitp = log_emit[cur][obs_idx[t]];
            if (emitp == NEG_INF) { dp_curr[cur] = NEG_INF; backptr[t][cur] = -1; continue; }
            double best = NEG_INF;
            int best_prev = -1;
            for (int prev=0; prev<S; ++prev) {
                if (dp_prev[prev] == NEG_INF) continue;
                double cand = dp_prev[prev] + log_trans[prev][cur];
                if (cand == NEG_INF) continue;
                if (cand > best) { best = cand; best_prev = prev; }
            }
            if (best_prev != -1) {
                dp_curr[cur] = best + emitp;
                backptr[t][cur] = best_prev;
            } else {
                dp_curr[cur] = NEG_INF;
                backptr[t][cur] = -1;
            }
        }
        for (int i=0;i<S;++i) dp_prev[i] = dp_curr[i];
    }

    double best_score = NEG_INF;
    int best_state = -1;
    for (int s=0;s<S;++s) {
        if (dp_prev[s] > best_score) { best_score = dp_prev[s]; best_state = s; }
    }
    if (best_state == -1) { cerr<<"No valid path found\n"; return 1; }

    int *path = new int[T];
    int cur = best_state;
    for (int t=T-1; t>=0; --t) {
        path[t] = cur;
        cur = backptr[t][cur];
    }

    cout<<fixed<<setprecision(6);
    cout<<"Most likely state sequence:\n";
    for (int t=0;t<T;++t) {
        cout << state_names[path[t]] << (t+1<T ? " " : "\n");
    }
    if (best_score <= -1e200) {
        cout << "Log-probability (approx): -inf\n";
    } else {
        cout << "Log-probability: " << best_score << "\n";
        double prob = (best_score < -700.0) ? 0.0 : exp(best_score);
        cout << "Probability (approx): " << prob << "\n";
    }

    return 0;
}
