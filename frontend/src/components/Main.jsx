import React from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route } from "react-router-dom";
import Header from './Header.jsx';
import Footer from './Footer.jsx';
import upload_contigs from './UploadContigs.jsx';
import Profile_view from './ProfileView.jsx';
import Upload_profile from './UploadProfile.jsx';
import Dendrogram_view from './DendrogramView.jsx';

class Main extends React.Component {

    render() {
        return(
            <Router>
                <div>
                    <Header />
                    <Route path="/" exact component={upload_contigs} />
                    <Route path="/profile_view" component={Profile_view} />
                    <Route path="/upload_profile" component={Upload_profile} />
                    <Route path="/dendrogram_view" component={Dendrogram_view} />
                    <Footer />
                </div>
            </Router>
        );
    }
}

ReactDOM.render(<Main /> ,document.getElementById('root'));
