import React from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route, HashRouter } from "react-router-dom";
import Header from './Header.jsx';
import Navigation from './Navigation.jsx'
import Footer from '../components/Footer.jsx';
import About from '../components/About.jsx';
import upload_contigs from './UploadContigs.jsx';
import Profile_view from './ProfileView.jsx';
import Tracking from './Tracking.jsx';
import Tracking_result from './TrackingResult.jsx';
import Upload_profile from './UploadProfile.jsx';
import Dendrogram_view from './DendrogramView.jsx';


class Main extends React.Component {

    constructor(props) {

        super(props);
        
        history.pushState(null, null, location.href);
        window.onpopstate = function(){
            history.go(1);
        };
    }

    render() {
        return(
            <Router>
                <div>
                    <Header />
                    <Navigation />
                    <Route path="/cgMLST/non-release/" exact component={upload_contigs} />
                    <Route path="/cgMLST/non-release/about" component={About} />
                    <Route path="/cgMLST/non-release/profiling" exact component={upload_contigs} />
                    <Route path="/cgMLST/non-release/profile_result" component={Profile_view} />
                    <Route path="/cgMLST/non-release/tracking" component={Tracking} />
                    <Route path="/cgMLST/non-release/tracking_result" component={Tracking_result} />
                    <Route path="/cgMLST/non-release/clustering" component={Upload_profile} />
                    <Route path="/cgMLST/non-release/clustering_result" component={Dendrogram_view} />
                    <Footer />
                </div>
            </Router>
        );
    }
}

ReactDOM.render(<Main />, document.getElementById('root'));
