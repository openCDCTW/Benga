import React from 'react';
import { Link } from 'react-router-dom';
import { withStyles } from '@material-ui/core/styles';
import List from '@material-ui/core/List';
import { ListItem, ListItemIcon, ListItemText } from '@material-ui/core/';
import InfoIcon from '@material-ui/icons/Info';
import FingerprintIcon from '@material-ui/icons/Fingerprint';
import TrackingIcon from '@material-ui/icons/GpsFixed';
import SearchIcon from '@material-ui/icons/Search';
import TreeIcon from '@material-ui/icons/AccountTree';
import FolderIcon from '@material-ui/icons/Folder';

const styles = theme => ({
  root: {
    marginTop: '20px',
  },
});

class Lists extends React.Component {

    constructor(props) {
        super(props);
    }

    render() {
        const { classes } = this.props;

        if(this.props.species == 'Vibrio_cholerae'){
            return(
                <div className={classes.root}>
                    <List>
                        <ListItem button component={Link} to="/cgMLST/about">
                        <ListItemIcon><InfoIcon /></ListItemIcon>
                        <ListItemText>About</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/profiling">
                        <ListItemIcon><FingerprintIcon /></ListItemIcon>
                        <ListItemText>cgMLST Profiling</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/tracking">
                        <ListItemIcon><TrackingIcon /></ListItemIcon>
                        <ListItemText>Strain Tracking</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/trackingSearch" disabled>
                        <ListItemIcon><SearchIcon /></ListItemIcon>
                        <ListItemText>Database search</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/clustering">
                        <ListItemIcon><TreeIcon /></ListItemIcon>
                        <ListItemText>Clustering</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/queryByID">
                        <ListItemIcon><FolderIcon /></ListItemIcon>
                        <ListItemText>Query result by ID</ListItemText>
                        </ListItem>
                    </List>
                </div>
            );
        }else if(this.props.species == 'Salmonella_enterica'){
            return(
                <div className={classes.root}>
                    <List>
                        <ListItem button component={Link} to="/cgMLST/profiling">
                        <ListItemIcon><FingerprintIcon /></ListItemIcon>
                        <ListItemText>cgMLST Profiling</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/clustering">
                        <ListItemIcon><TreeIcon /></ListItemIcon>
                        <ListItemText>Clustering</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/queryByID">
                        <ListItemIcon><FolderIcon /></ListItemIcon>
                        <ListItemText>Query result by ID</ListItemText>
                        </ListItem>
                    </List>
                </div>
            );
        }else if(this.props.species == 'Neisseria_meningitidis'){
            return(
                <div className={classes.root}>
                    <List>
                        <ListItem button component={Link} to="/cgMLST/about">
                        <ListItemIcon><InfoIcon /></ListItemIcon>
                        <ListItemText>About</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/profiling">
                        <ListItemIcon><FingerprintIcon /></ListItemIcon>
                        <ListItemText>cgMLST Profiling</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/tracking">
                        <ListItemIcon><TrackingIcon /></ListItemIcon>
                        <ListItemText>Strain Tracking</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/trackingSearch" disabled>
                        <ListItemIcon><SearchIcon /></ListItemIcon>
                        <ListItemText>Database search</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/clustering">
                        <ListItemIcon><TreeIcon /></ListItemIcon>
                        <ListItemText>Clustering</ListItemText>
                        </ListItem>
                        <ListItem button component={Link} to="/cgMLST/queryByID">
                        <ListItemIcon><FolderIcon /></ListItemIcon>
                        <ListItemText>Query result by ID</ListItemText>
                        </ListItem>
                    </List>
                </div>
            );
        }
    }
}

export default withStyles(styles)(Lists);