import React from 'react';
import ReactDOM from 'react-dom';
import PropTypes from 'prop-types';
import { withStyles } from '@material-ui/core/styles';
import Input from '@material-ui/core/Input';
import InputBase from '@material-ui/core/InputBase';
import InputLabel from '@material-ui/core/InputLabel';
import TextField from '@material-ui/core/TextField';
import FormControl from '@material-ui/core/FormControl';
import SearchIcon from '@material-ui/icons/Search';
import Button from '@material-ui/core/Button';
import { Link } from 'react-router-dom';

const styles = theme => ({
  container: {
    display: 'flex',
    flexWrap: 'wrap',
  },
  margin: {
    margin: theme.spacing.unit,
  },
  bootstrapRoot: {
    'label + &': {
      marginTop: theme.spacing.unit * 3,
    },
  },
  bootstrapInput: {
    borderRadius: 4,
    backgroundColor: theme.palette.common.white,
    border: '1px solid #ced4da',
    fontSize: 16,
    padding: '10px 12px',
    transition: theme.transitions.create(['border-color', 'box-shadow']),
    // Use the system font instead of the default Roboto font.
    '&:focus': {
      borderColor: '#80bdff',
      boxShadow: '0 0 0 0.2rem rgba(0,123,255,.25)',
    },
  },
  bootstrapFormLabel: {
    fontSize: 18,
  },
});

class InputBatchID extends React.Component {

  handleChange = ID => event => {
    window.queryID = event.target.value
    }

  render() {
    const { classes } = this.props;

    return (
      <div className={classes.container}>
        <FormControl className={classes.margin}>
          <InputLabel shrink className={classes.bootstrapFormLabel}>
            Batch ID
        </InputLabel>
          <InputBase
            placeholder="Input your batch ID here"
            classes={{
              root: classes.bootstrapRoot,
              input: classes.bootstrapInput,
            }}
            onChange={this.handleChange('ID')}
          />
        </FormControl>
        <Link to="/query_data" style={{ textDecoration:'none',display:'flex' }}>
          <Button variant="outlined" color="default" size="medium">
              <SearchIcon />
          </Button>
        </Link>
      </div>
    );
  }
}

InputBatchID.propTypes = {
  classes: PropTypes.object.isRequired,
};

export default withStyles(styles)(InputBatchID);
